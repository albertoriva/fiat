#!/usr/bin/env python

import sys
import curses
import importlib

if sys.version_info.major == 2:
    cw = importlib.import_module("curses.wrapper")
    wrapper = cw.wrapper
else:
    wrapper = curses.wrapper

__doc__       = "Text-mode browser for DNA sequences"
__author__    = "Alberto Riva"
__version__   = "0.9"
__copyright__ = "(c) 2020, A. Riva, University of Florida"


BMATCHES = ["KG", "KT", "GK", "TK",               # G/T
            "YC", "YT", "CY", "TY",               # C/T
            "SC", "SG", "CS", "GS",               # C/G
            "WA", "WT", "AW", "TW",               # A/T
            "RA", "RG", "AR", "GR",               # A/G
            "MA", "MC", "AM", "CM",               # A/C
            "BC", "BG", "BT", "CB", "GB", "TB",   # C/G/T
            "DA", "DG", "DT", "AD", "GD", "TD",   # A/G/T
            "HA", "HC", "HT", "AH", "CH", "TH",   # A/C/T
            "VA", "VC", "VG", "AV", "CV", "GV"]   # A/C/G

VALIDBASES = "ACGTKYSWRDMHVNX"

BASECOMPLEMENTS = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                   'K': 'M', 'Y': 'R', 'S': 'S',
                   'W': 'W', 'R': 'Y', 'M': 'K',
                   'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B'}

ORIENTATIONS = {'f': 'fwd',
                'r': 'rev',
                'c': 'comp',
                'd': 'revcomp'}

REGCOLORS = {'r': 3, 'g': 4, 'b': 5, 'm': 6, 'c': 7}

def bmatch(b1, b2):
    if b1 == b2:
        return True
    if b1 == 'N' or b2 == 'N':
        return True
    if b1+b2 in BMATCHES:
        return True
    return False

def seqmatch(shortseq, sl, longseq, ll, p, mm=0):
    """Returns True if `shortseq' of length `sl' matches sequence `longseq' of length `ll'
starting at position p, with at most mm mismatches."""
    maxp = len(longseq)
    nm = 0
    for i in range(sl):
        if not bmatch(shortseq[i], longseq[p]):
            nm += 1
            if nm > mm:
                return False
        p += 1
        if p > ll:
            return False
    return True

def validateSeq(seq):
    for b in seq:
        if b not in VALIDBASES:
            return False
    return True

def validSearchMode(sm):
    result = ""
    for m in "frcd":
        if m in sm:
            result += m
    return result or "f"

def decodeOrientation(o):
    if o in ORIENTATIONS:
        return ORIENTATIONS[o]
    else:
        return o

def identity(seq):
    return seq
    
def reverseSeq(seq):
    return seq[::-1]

def complementSeq(seq):
    return "".join([ BASECOMPLEMENTS[b] for b in seq ])

def reverseComplementSeq(seq):
    return "".join([ BASECOMPLEMENTS[b] for b in seq[::-1] ])

ORIENTFUNC = {'f': identity,
              'r': reverseSeq,
              'c': complementSeq,
              'd': reverseComplementSeq}

def decodeColor(s):
    s = s.lower()
    if s in REGCOLORS:
        return REGCOLORS[s]
    else:
        return REGCOLORS['r']

def parseCoords(c):
    try:
        if "-" in c:
            [a, b] = c.split("-")
            return [int(a), int(b)]
        elif "+" in c:
            [a, b] = c.split("+")
            return [int(a), int(a) + int(b)]
        else:
            return [False, False]
    except ValueError:
        return [False, False]
    
class Region(object):
    seq = ""
    regname = ""
    orientation = "f"           # Or "r", "c", "d"
    color = 3
    start = -1
    end = -1

    def __init__(self, start, end, seq, prefix='r', color=3, name=None):
        self.seq = seq
        self.start = start
        self.end = end
        if name:
            self.regname = name
        else:
            self.regname = "{}_{}_{}".format(prefix, start, end)
        self.color = color
    
class Sequence(object):
    fastafile = ""
    seq = ""
    seqname = ""
    seqlen = 0
    label = ""
    regions = []
    regidx = 0
    hits = []
    hitidx = -1
    searchmode = "f"
    mismatches = 0
    stats = {}
    status = None
    _focus = None               # Hit or region being displayed
    _focustype = ""             # 'h' or 'r'
    _nrows = 0
    _row = 0                    # First row being displayed
    _bottom = False             # True if we're displaying the last line in the sequence
    
    def __init__(self):
        self.regions = []
        self.hits = []

    def load(self, fastafile):
        self.fastafile = fastafile
        with open(fastafile, "r") as f:
            hdr = f.readline()
            if hdr[0] == '>':
                self.seqname = hdr.rstrip("\r\n")[1:]
            else:
                sys.stderr.write("Malformed FASTA file - cannot load.\n")
                return False
            for line in f:
                self.seq += line.rstrip().upper()
        self.seqlen = len(self.seq)
        self._nrows = 1 + self.seqlen / 60
        self.calcStats()
        sys.stderr.write("Sequence {} loaded, {}bp.\n".format(self.seqname, self.seqlen))
        return True

    def finalize(self, label):
        self.label = label
        self.seqlen = len(self.seq)
        self._nrows = 1 + self.seqlen / 60
        self.calcStats()
    
    def calcStats(self):
        for b in "ACGT":
            self.stats[b] = self.seq.count(b)
    
    def validRow(self, r):
        if r < 0:
            return 0
        elif r > self._nrows:
            return self._nrows
        else:
            return r

    def findMatches(self, target):
        tl = len(target)
        targets = []
        for m in self.searchmode:
            func = ORIENTFUNC[m]
            targets.append((func(target), m))
        self.hits = []
        for p in range(0, self.seqlen - tl):
            for tg in targets:
                if seqmatch(tg[0], tl, self.seq, self.seqlen, p, mm=self.mismatches):
                    reg = Region(p+1, p+tl, self.seq[p:p+tl], prefix='h', color=2)
                    reg.orientation = tg[1]
                    self.hits.append(reg)

    # Regions and hits

    def getOverlapping(self, p, regions):
        oh = []
        for h in regions:
            if h.start > p:
                break
            if h.start <= p <= h.end:
                oh.append(h)
        return oh

    def getAllOverlapping(self, start, end, regions):
        oh = []
        for h in regions:
            if h.start > end:
                break
            if start <= h.start < end or start <= h.end < end or (h.start <= start and h.end >= end):
                oh.append(h)
        return oh
    
    def getHitsAt(self, p):
        """Returns a list of all hits overlapping position `p'."""
        return self.getOverlapping(p, self.hits)

    def getRegionsAt(self, p):
        """Returns a list of all regions overlapping position `p'."""
        return self.getOverlapping(p, self.regions)

    def getHitsRange(self, start, end):
        return self.getAllOverlapping(start, end, self.hits)

    def getRegionsRange(self, start, end):
        return self.getAllOverlapping(start, end, self.regions)

    # Interface methods

    #+IGNORE
    def getSequence(self, p):
        blocks = []
        hits = self.getHitsRange(p, p+60)
        if hits:
            starts = sorted([r.start for r in hits])
            ends   = sorted([r.end+1 for r in hits])
            openh  = len([e for e in starts if e < p])
            bnames = [h.regname for h in hits]
            #print (starts, ends, openh)
        else:
            starts = []
            ends   = []
            openh  = 0
            bnames = []
        bstart = p
        mode   = "H" if openh > 0 else "N"             # H = inside hits, N = no hits
        for i in range(p, p+61):
            openb = openh
            if i in starts:
                openh += 1
            if i in ends:
                openh += -1
            if mode == "H":
                if openh == 0:
                    if i > bstart:
                        blocks.append((bstart, i-1, 1))
                    bstart = i
                    mode = "N"
            else:
                if openh > 0:
                    if i > bstart:
                        blocks.append((bstart, i-1, 0))
                    bstart = i
                    mode = "H"
        if openh:
            blocks.append((bstart, i-1, 1))
        else:
            blocks.append((bstart, i-1, 0))
        return (blocks, bnames)

    def getBlocks(self, p):
        blocks = []
        bnames = []
        buf = [0]*60
        pe = p+60
        regions = self.getRegionsRange(p, pe)
        hits = self.getHitsRange(p, pe)
        for reg in regions:
            for i in range(reg.start, reg.end+1):
                ii = i-p
                if 0 <= ii < 60:
                    buf[ii] = reg.color
        for hit in hits:
            for i in range(hit.start, hit.end+1):
                ii = i - p
                if 0 <= ii < 60:
                    buf[ii] = 1
        for hit in hits:
            bnames.append((hit.regname, hit.color))
        for reg in regions:
            bnames.append((reg.regname, reg.color))
        bstart = p
        bend   = p + 1
        col = buf[0]
        for i in range(1, 60):
            if buf[i] != col:
                blocks.append((bstart, bend-1, col))
                bstart = bend
                col = buf[i]
            bend += 1
        blocks.append((bstart, pe - 1, col))
        return (blocks, bnames)

    def display(self, win):
        self._bottom = False
        (h, w) = win.getmaxyx()
        maxrow = h - 2
        r = self._row           # First row to display
        p = r * 60 + 1          # First position to display (1-based)
        start = p
        win.move(0, 0)
        win.erase()
        win.addstr(0, 0, "                    1         2         3         4         5         6")
        win.addstr(1, 0, "           123456789012345678901234567890123456789012345678901234567890")
        ypos = 3
        while True:
            win.addstr(ypos, 0, "{:10} ".format(p))
            (blocks, bnames) = self.getBlocks(p) # self.getSequence(p)
            for bl in blocks:
                win.addstr(self.seq[bl[0]-1:bl[1]], curses.color_pair(bl[2]))
            win.addstr(" ")
#            win.addstr(" {}".format(blocks))
#            win.addstr(" {}".format(bnames))
            for bn in bnames:
                win.addstr(" " + bn[0], curses.color_pair(bn[1]))
            ypos += 1
            r += 1
            p += 60
            if p >= self.seqlen:
                self._bottom = True
                break
            if ypos == maxrow:
                break
        statusline = ">{} | Length: {}bp | Range: {}-{} | Regions: {} | Hits: {}".format(self.seqname, self.seqlen, start, p-1, len(self.regions), len(self.hits))
        statusline2 = self.label
        nsp = w - len(statusline) - len(statusline2)
        win.addstr(maxrow, 0, statusline + " "*nsp + statusline2, curses.A_REVERSE)
        if not self.status:
            self.status = "Press '?' for help."
        win.addstr(maxrow+1, 0, self.status)
        self.status = None
        win.refresh()

    def right(self):
        self._row = min(self._row + 10, self._nrows)

    def left(self):
        self._row = max(0, self._row - 10)

    def down(self):
        if not self._bottom:
            self._row += 1

    def up(self):
        if self._row > 0:
            self._row += -1

    def pageup(self, win):
        (h, w) = win.getmaxyx()
        h = h -5
        if self._row >= h:
            self._row -= h

    def pagedown(self, win):
        (h, w) = win.getmaxyx()
        h = h - 5
        if self._row + h < self._nrows:
            self._row += h

    def top(self):
        self._row = 0

    def bottom(self):
        self._row = self._nrows - 1

    def askInt(self, win, prompt):
        (h, w) = win.getmaxyx()
        win.move(h-1, 0)
        win.clrtoeol()
        win.addstr(prompt, curses.A_BOLD)
        try:
            curses.echo()
            curses.curs_set(1)
            r = win.getstr()
        finally:
            curses.curs_set(0)
            curses.noecho()
        try:
            return int(r)
        except ValueError:
            return None

    def askString(self, win, prompt):
        (h, w) = win.getmaxyx()
        win.move(h-1, 0)
        win.clrtoeol()
        win.addstr(prompt, curses.A_BOLD)
        try:
            curses.echo()
            curses.curs_set(1)
            r = win.getstr()
        finally:
            curses.curs_set(0)
            curses.noecho()
        return r

    def askColor(self, win):
        col = self.askString(win, "Color (r,g,b,c,m): ")
        return decodeColor(col)
    
    def showMessage(self, win, message, wait=True):
        (h, w) = win.getmaxyx()
        win.move(h-1, 0)
        win.clrtoeol()
        win.addstr(message, curses.A_BOLD)
        win.refresh()
        if wait:
            return win.getch()
        else:
            return False
        
    def askPosition(self, win):
        r = self.askInt(win, "Go to position: ")
        if not r:
            return
        self.goto(r)
        
    def goto(self, r, delta=0):
        if r <= 0 or r > self.seqlen:
            return
        self._row = self.validRow(r / 60 + delta)
        
    def doSearch(self, win):
        target = self.askString(win, "Search: ")
        if not target:
            return
        if target == "-":
            self.hits = []
            self.hitidx = -1
            return
        target = target.upper()
        if not validateSeq(target):
            self.showMessage(win, "Error: search sequence contains invalid nucleotides (allowed: {})".format(VALIDBASES))
            return
        self.findMatches(target)
        self.status = "Search: {} hits found.".format(len(self.hits))
        
    def nextHit(self, win):
        nh = len(self.hits)
        if nh == 0:
            return
        self.hitidx += 1
        if self.hitidx == nh:
            self.hitidx = 0
        return self.showHit(win, nh)
        
    def prevHit(self, win):
        nh = len(self.hits)
        if nh == 0:
            return
        self.hitidx += -1
        if self.hitidx == -2:
            self.hitidx = 0
        elif self.hitidx == -1:
            self.hitidx = nh - 1
        return self.showHit(win, nh)

    def showHit(self, win, nh):
        hit = self.hits[self.hitidx]
        self.goto(hit.start, -1)
        self.status = "* Hit {}/{} | Position: {}-{} | Orientation: {} | (a)dd as region | (A)dd all as regions".format(self.hitidx+1, nh, hit.start, hit.end, decodeOrientation(hit.orientation))
        return hit

    def nextRegion(self, win):
        nh = len(self.regions)
        if nh == 0:
            return
        self.regidx += 1
        if self.regidx == nh:
            self.regidx = 0
        return self.showRegion(win, nh)
        
    def prevRegion(self, win):
        nh = len(self.regions)
        if nh == 0:
            return
        self.regidx += -1
        if self.regidx == -2:
            self.regidx = 0
        elif self.regidx == -1:
            self.regidx = nh - 1
        return self.showRegion(win, nh)

    def showRegion(self, win, nh):
        hit = self.regions[self.regidx]
        self.goto(hit.start, -1)
        self.status = "* Region {}/{} | Position: {}-{} | Orientation: {} | (r)ename | (c)olor | (d)elete".format(self.regidx+1, nh, hit.start, hit.end, decodeOrientation(hit.orientation))
        return hit

    def renameRegion(self, win):
        newname = self.askString(win, "New name: ")
        if newname:
            self._focus.regname = newname

    def recolorRegion(self, win):
        newcolor = self.askColor(win)
        self._focus.color = newcolor
            
    def addRegion(self, win):
        coords = self.askString(win, "Enter position (start-end or start+length): ")
        (start, end) = parseCoords(coords)
        if start and end:
            regname = "r_{}_{}".format(start, end)
            newname = self.askString(win, "Region name [{}]: ".format(regname))
            if newname:
                regname = newname
            color = self.askColor(win)
            reg = Region(start, end, self.seq[start:end+1], color, name=regname)
            self.regions.append(reg)

    def deleteRegion(self, win):
        self.regions.remove(self._focus)
        self._focus = None
            
    def deleteAllRegions(self, win):
        ans = self.askString(win, "Delete all regions (yes/no)? ")
        if ans == "yes":
            self.regions = []
            self._focus = None
            
    def hitToRegion(self, win):
        color = self.askColor(win)
        self.hits.remove(self._focus)
        self.regions.append(self._focus)
        self._focus.color = color
        self.regions.sort(key=lambda r: r.start)

    def allHitsToRegions(self, win):
        color = self.askString(win, "Color (r,g,b,c,m): ")
        col = decodeColor(color)
        for hit in self.hits:
            hit.color = col
            self.regions.append(hit)
        self.hits = []
        self.regions.sort(key=lambda r: r.start)
    
    def showStats(self, win):
        msg = "Length: {}bp, A: {} ({:.1f}%), C: {} ({:.1f}%), G: {} ({:.1f}%), T: {} ({:.1f}%), GC%: {:.1f}".format(
            self.seqlen, self.stats['A'], 100.0 * self.stats['A'] / self.seqlen,
            self.stats['C'], 100.0 * self.stats['C'] / self.seqlen,
            self.stats['G'], 100.0 * self.stats['G'] / self.seqlen,
            self.stats['T'], 100.0 * self.stats['T'] / self.seqlen,
            100.0 * (self.stats['C'] + self.stats['G']) / self.seqlen)
        self.showMessage(win, msg)

    def showOptions(self, win):
        while True:
            a = self.showMessage(win, "Options: [m]ismatches = {} | [s]earch mode = {}".format(self.mismatches, self.searchmode))
            if a == ord('m'):
                mm = self.askInt(win, "Max mismatches in search: ")
                if mm:
                    self.mismatches = mm
            elif a == ord('s'):
                sm = self.askString(win, "Search mode (one or more of f,r,c,d): ")
                if sm:
                    sm = validSearchMode(sm)
                    self.searchmode = sm
            else:
                return
                
class Driver(object):
    fastas = []
    nfastas = 0
    cmd = None
    
    def __init__(self, args):
        self.fastas = []
        self.nfastas = 0
        self.parseArgs(args)
        
    def parseArgs(self, args):
        for a in args:
            self.addSequences(a)
        self.nfastas = len(self.fastas)
        idx = 1
        for fa in self.fastas:
            fa.finalize("[{}/{}]".format(idx, self.nfastas))
            idx += 1

    def addSequences(self, filename):
        with open(filename, "r") as f:
            for line in f:
                if line[0] == '>':
                    S = Sequence()
                    S.seqname = line.rstrip("\r\n")[1:]
                    self.fastas.append(S)
                else:
                    S.seq += line.rstrip().upper()
        
    def displayAll(self, win):
        curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_YELLOW)  # search hit in seq
        curses.init_pair(2, curses.COLOR_YELLOW, curses.COLOR_BLACK)  # search hit name
        curses.init_pair(3, curses.COLOR_RED, curses.COLOR_BLACK)     # r
        curses.init_pair(4, curses.COLOR_GREEN, curses.COLOR_BLACK)   # g
        curses.init_pair(5, curses.COLOR_BLUE, curses.COLOR_BLACK)    # b
        curses.init_pair(6, curses.COLOR_MAGENTA, curses.COLOR_BLACK) # m
        curses.init_pair(7, curses.COLOR_CYAN, curses.COLOR_BLACK)    # c
        
        fidx = 0
        while True:
            self.run(win, self.fastas[fidx])
            if self.cmd == 'quit':
                return
            elif self.cmd == 'next':
                fidx += 1
                if fidx == self.nfastas:
                    fidx = 0
            elif self.cmd == 'prev':
                fidx += -1
                if fidx < 0:
                    fidx = self.nfastas - 1
            self.cmd = None
            
    def run(self, win, fasta):
        curses.curs_set(0)
        while True:
            fasta.display(win)
            a = win.getch()
            if a in [113, 81]:      # Quit (q, Q)
                self.cmd = 'quit'
                break
            elif a == ord('n'):
                self.cmd = 'next'
                break
            elif a == ord('p'):
                self.cmd = 'prev'
                break
            elif a == ord('?'):
                self.showHelp(win)
            elif a == curses.KEY_RIGHT:
                fasta.right()
            elif a == curses.KEY_LEFT:
                fasta.left()
            elif a == curses.KEY_UP:
                fasta.up()
            elif a in [curses.KEY_DOWN, ord('\n')]:
                fasta.down()
            elif a == curses.KEY_HOME:
                fasta.top()
            elif a == curses.KEY_END:
                fasta.bottom()
            elif a == curses.KEY_PPAGE:
                fasta.pageup(win)
            elif a in [curses.KEY_NPAGE, 32]:
                fasta.pagedown(win)
            elif a == ord('g'):
                fasta.askPosition(win)
            elif a == ord('s'):
                fasta.showStats(win)
            elif a == ord('/'):
                fasta.doSearch(win)
            elif a == ord('o'):
                fasta.showOptions(win)
            elif a == ord(','):
                fasta._focustype = 'h'
                fasta._focus = fasta.prevHit(win)
            elif a == ord('.'):
                fasta._focus = fasta.nextHit(win)
                fasta._focustype = 'h'
            elif a == ord('<'):
                fasta._focus = fasta.prevRegion(win)
                fasta._focustype = 'r'
            elif a == ord('>'):
                fasta._focus = fasta.nextRegion(win)
                fasta._focustype = 'r'
            elif a == ord('a'):
                if fasta._focus and fasta._focustype == 'h':
                    fasta.hitToRegion(win)
                else:
                    fasta.addRegion(win)
            elif a == ord('A'):
                if fasta._focus and fasta._focustype == 'h':
                    fasta.allHitsToRegions(win)
            elif a == ord('r'):
                if fasta._focus and fasta._focustype == 'r':
                    fasta.renameRegion(win)
            elif a == ord('c'):
                if fasta._focus and fasta._focustype == 'r':
                    fasta.recolorRegion(win)
            elif a == ord('d'):
                if fasta._focus and fasta._focustype == 'r':
                    fasta.deleteRegion(win)
            elif a == ord('D'):
                fasta.deleteAllRegions(win)
                    
    def showHelpPage(self, win, text, footer=None):
        (h, w) = win.getmaxyx()
        win.move(0, 0)
        win.erase()
        win.addstr(1, 1, "FIAT - FASTA In A Terminal", curses.color_pair(4) + curses.A_BOLD)
        win.addstr(3, 0, text)
        if footer:
            win.addstr(h-2, 0, footer)
        win.addstr(h-1, 0, "(c) 2020, A. Riva, University of Florida", curses.A_DIM)
        return win.getch()
    
    def showHelp(self, win):
        a = self.showHelp1(win)
        if a == ord('q'):
            return
        a = self.showHelp2(win)
        if a == ord('q'):
            return
        a = self.showHelp3(win)

    def showHelp1(self, win):
        return self.showHelpPage(win, """
        FIAT displays one or more FASTA sequences in a terminal window. The sequence
        is displayed across the whole height of the terminal except for the top three
        lines (used for coordinates) and the bottom two. The last-but-one row (in
        reverse) displays the sequence name and length, the currently visible coordinate
        range, the number of regions and the number of hits. At the right end of the
        line, it shows the number of the current sequence and the total number of sequences.
        The last line (the `message' line) is used to display transient information or
        to get use input.

 Basic commands:

        Up, Down    - scroll one line backwards / forward
        Left, Right - scroll 10 lines backwards / forward
        PgUp, PgDn  - scroll one page backwards / forward
        Space       - scroll one page forward
        Home, End   - go to beginning / end of sequence
        g           - jump to specified position
        s           - display sequence statistics
        o           - display options
        n, p        - switch to next / previous sequence
        /           - search for patterns
        a           - add a region
        D           - delete all regions
        q           - quit

 When displaying options:

        m           - set maximum number of mismatches in search
        s           - set search mode - any subset of the following:
                        f = forward, r = reverse, c = complement, d = reverse complement

        """, "Press q to exit help, any other key for next page (Regions).")

    def showHelp2(self, win):
        return self.showHelpPage(win, """Regions

        Regions are arbitrary subsequences characterized by a start and end position,
        a name, and a color. Regions can be defined manually (using the `a' key) or
        as a result of a search. Use < and > to focus the previous / next region
        respectively. When a region is focused, the message line starts with `* Region'
        and shows the number of the current region and its coordinates.

        <, >         - jump to previous / next region
        r            - rename the currently focused region
        c            - change color of the currently focused region
        d            - delete the currently focused region

        The following five colors can be used for regions, identified by their initial:
          (r)ed, (g)reen, (b)lue, (m)agenta, (c)yan

        """, "Press q to exit help, any other key for next page (Searching).")

    def showHelp3(self, win):
        return self.showHelpPage(win, """Searching

        The search command finds all occurrences of a specified pattern (hits). Hits are
        displayed in black on yellow in the sequence, and their names are displayed in
        yellow. The name of a search hit is automatically generated and cannot be changed. 

        The pattern can be specified using the full IUPAC nucleotide alphabet, shown here:

        A      K = G/T       B = C/G/T
        C      Y = C/T       D = A/G/T
        G      S = C/G       H = A/C/T
        T      W = A/T       V = A/C/G
               R = A/G
               M = A/C       N = A/C/G/T

        For example, the pattern CAST will match CACT or CAGT.

        Matching is also affected by the maximum number of mismatches allowed (that can be
        set using the `o' command followed by `m') and by the search mode (set with `o'
        followed by `s'), which can be any subset of the following:

          f = forward, r = reverse, c = complement, d = reverse complement
        
        For example, if the search mode is fd, GATTA will match both GATTA and TAATC.

        Use , and . to focus the previous / next hit respectively. When a hit is focused, 
        the message line starts with `* Hit' and shows the number of the current hit, its 
        coordinates, and its orientation.

        .            - jump to hte previous hit
        ,            - jump to the next hit
        a            - save this hit as a region (prompts for name and color)
        A            - save all hits as regions (prompts for color)

        To clear the list of hits, enter `-' as the search string.

        """, "Press any key to exit help.")
# Main

def usage():
    sys.stdout.write("""fiat.py - FASTA In A Terminal

Usage: fiat.py fastafiles...

Displays one or more sequences from FASTA files in the terminal. Allows
scrolling through the sequence, searching for subsequences, defining
arbitrary regions identified by a name and a color.

Press the '?' key when displaying a sequence for detailed instructions.

{}

""".format(__copyright__))
    
if __name__ == "__main__":
    args = sys.argv[1:]
    if "-h" in args or "--help" in args:
        usage()
    else:
        D = Driver(args)
        if D.nfastas > 0:
            wrapper(D.displayAll)
        else:
            usage()     
