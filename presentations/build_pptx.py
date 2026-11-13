#!/usr/bin/env python3
"""Build presentations/all-the-data-deck.pptx from presentations/all-the-data-deck.html.

The HTML deck is the source. Every <section> becomes one slide, in order: claim title, kicker,
bullets, boxes, tables, diagrams, source line and speaker notes (<aside class="notes"> goes to the
notes pane). The slides of the section "Update, 12 November 2026" are written in UPDATE_HTML below,
in the same markup, and inserted before the appendix.

Run from the repository root:

  # 1. Diagrams: render diagrams/d*.svg to diagrams/png/d*.png with the installed Google Chrome.
  #    The SVGs carry HTML labels (foreignObject), which rsvg-convert and Inkscape do not draw.
  uv run --no-project --with playwright==1.48.0 python presentations/build_pptx.py --diagrams

  # 2. The deck (reads the PNGs from step 1).
  uv run --no-project --with python-pptx==1.0.2 --with beautifulsoup4==4.12.3 \\
      python presentations/build_pptx.py

Text heights are estimated with the Arial metrics found on the host (macOS or msttcorefonts); each
slide takes the largest font sizes that fit (body 18 to 16 pt, tables 18 to 12 pt) and gives the
diagram the height that is left.
"""

from __future__ import annotations

import argparse
import datetime as dt
import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

HERE = Path(__file__).resolve().parent
HTML = HERE / "all-the-data-deck.html"
OUT = HERE / "all-the-data-deck.pptx"
DIAG = HERE / "diagrams"
PNG = DIAG / "png"
CHROME = "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"

# Palette of the HTML deck (:root variables).
INK, MUTED, RULE = "14213D", "55607A", "D6DBE5"
BLUE, BLUE_BG = "2F5D9B", "E8EEF7"
GREEN, GREEN_BG = "3C7D40", "EEF6EE"
AMBER, AMBER_BG = "C26A1B", "FDF1E3"
RED, RED_BG = "B23A48", "FBECEE"
PANEL, WHITE, PEACH, SLATE = "F7F8FB", "FFFFFF", "F2B880", "55607A"
VARS = {"--ink": INK, "--muted": MUTED, "--rule": RULE, "--blue": BLUE, "--blue-bg": BLUE_BG,
        "--green": GREEN, "--green-bg": GREEN_BG, "--amber": AMBER, "--amber-bg": AMBER_BG,
        "--red": RED, "--red-bg": RED_BG, "--panel": PANEL}
BOX = {None: (PANEL, RULE), "blue": (BLUE_BG, BLUE), "green": (GREEN_BG, GREEN),
       "amber": (AMBER_BG, AMBER), "red": (RED_BG, RED)}
TAG = {"pod": GREEN, "edge": AMBER, "off": BLUE, "lic": RED}

FONT, MONO = "Arial", "Courier New"
SLIDE_W, SLIDE_H = 13.333, 7.5
X0, W0 = 0.55, 13.333 - 1.1          # content left edge and width (inches)
TITLE_PT, KICKER_PT, SMALL_PT, SRC_PT = 24, 14, 14, 12
LH = 1.2                             # line height as a multiple of the font size
SPACE_AFTER = 3                      # points after each paragraph
MIN_IMG = 3.0                      # inches: smallest diagram height before fonts shrink


# --------------------------------------------------------------------------- diagrams

def render_diagrams(scale: int = 3) -> None:
    """Screenshot each Mermaid SVG at `scale` times its own size, on white."""
    from playwright.sync_api import sync_playwright

    PNG.mkdir(exist_ok=True)
    with sync_playwright() as p:
        browser = p.chromium.launch(executable_path=CHROME)
        page = browser.new_context(device_scale_factor=scale).new_page()
        for svg in sorted(DIAG.glob("d*.svg")):
            page.goto(svg.as_uri())
            w, h = page.evaluate(
                """() => { const s = document.documentElement, vb = s.viewBox.baseVal;
                   s.setAttribute('width', vb.width); s.setAttribute('height', vb.height);
                   s.style.maxWidth = 'none'; s.style.backgroundColor = '#ffffff';
                   return [Math.ceil(vb.width), Math.ceil(vb.height)]; }""")
            page.set_viewport_size({"width": w, "height": h})
            page.wait_for_timeout(200)
            out = PNG / f"{svg.stem}.png"
            page.screenshot(path=str(out), clip={"x": 0, "y": 0, "width": w, "height": h})
            print(f"{out.relative_to(HERE)}  {w * scale}x{h * scale}")
        browser.close()


# --------------------------------------------------------------------------- text metrics

_FONTS: dict = {}
_FONT_FILES = {(False, False): "Arial.ttf", (True, False): "Arial Bold.ttf",
               (False, True): "Arial Italic.ttf", (True, True): "Arial Bold Italic.ttf",
               "mono": "Courier New.ttf"}
_FONT_DIRS = [Path("/System/Library/Fonts/Supplemental"), Path("/Library/Fonts"),
              Path("/usr/share/fonts/truetype/msttcorefonts")]


def _font(key):
    if key not in _FONTS:
        from PIL import ImageFont
        _FONTS[key] = None
        for d in _FONT_DIRS:
            f = d / _FONT_FILES[key]
            if f.exists():
                _FONTS[key] = ImageFont.truetype(str(f), 100)
                break
    return _FONTS[key]


def text_w(s: str, size: float, bold=False, italic=False, code=False) -> float:
    """Width of `s` in points."""
    f = _font("mono" if code else (bold, italic))
    if f is None:
        return len(s) * size * (0.6 if code else 0.55)
    return f.getlength(s) * size / 100


# --------------------------------------------------------------------------- inline model

@dataclass
class Run:
    text: str
    bold: bool = False
    italic: bool = False
    code: bool = False
    color: str | None = None


@dataclass
class Para:
    runs: list
    kind: str = "p"          # p | h3 | li
    ordered: bool = False

    def text(self) -> str:
        return "".join(r.text for r in self.runs)


def css_color(value: str | None) -> str | None:
    if not value:
        return None
    value = value.strip()
    m = re.match(r"var\((--[\w-]+)\)", value)
    if m:
        return VARS.get(m.group(1))
    m = re.match(r"#([0-9a-fA-F]{6})", value)
    return m.group(1).upper() if m else None


def style_of(el) -> dict:
    out = {}
    for part in (el.get("style") or "").split(";"):
        if ":" in part:
            k, v = part.split(":", 1)
            out[k.strip()] = v.strip()
    return out


def finish(runs: list, kind="p", ordered=False) -> Para | None:
    merged: list[Run] = []
    for r in runs:
        if not r.text:
            continue
        if merged and merged[-1].text.endswith(" ") and r.text.startswith(" "):
            r = Run(r.text.lstrip(" "), r.bold, r.italic, r.code, r.color)
            if not r.text:
                continue
        prev = merged[-1] if merged else None
        if prev and (prev.bold, prev.italic, prev.code, prev.color) == (r.bold, r.italic, r.code, r.color):
            prev.text += r.text
        else:
            merged.append(Run(r.text, r.bold, r.italic, r.code, r.color))
    while merged and not merged[0].text.strip():
        merged.pop(0)
    while merged and not merged[-1].text.strip():
        merged.pop()
    if not merged:
        return None
    merged[0].text = merged[0].text.lstrip()
    merged[-1].text = merged[-1].text.rstrip()
    return Para(merged, kind, ordered)


def parse_nodes(nodes, kind="p") -> list[Para]:
    """Paragraphs of a run of HTML nodes: block elements break paragraphs, inline ones format runs."""
    from bs4 import Comment, NavigableString

    out: list[Para] = []
    cur: list[Run] = []

    def flush():
        nonlocal cur
        p = finish(cur, kind)
        cur = []
        if p:
            out.append(p)

    def one_para(node, fmt, k, ordered=False):
        sub = parse_nodes(node.children)
        runs = []
        for p in sub:
            for r in p.runs:
                runs.append(Run(r.text, r.bold or fmt.get("bold", False), r.italic, r.code,
                                r.color or fmt.get("color")))
            runs.append(Run(" "))
        return finish(runs, k, ordered)

    def walk(children, fmt):
        for c in children:
            if isinstance(c, Comment):
                continue
            if isinstance(c, NavigableString):
                t = re.sub(r"\s+", " ", str(c))
                if t:
                    cur.append(Run(t, **fmt))
                continue
            name, cls = c.name, (c.get("class") or [])
            if name == "aside":
                continue
            if name == "br":
                flush()
                continue
            if name in ("h1", "h2", "h3"):
                flush()
                p = one_para(c, fmt, "h3")
                if p:
                    out.append(p)
                continue
            if name in ("ul", "ol"):
                flush()
                for li in c.find_all("li", recursive=False):
                    p = one_para(li, fmt, "li", ordered=(name == "ol"))
                    if p:
                        out.append(p)
                continue
            if name in ("p", "div", "table", "tr", "td", "th", "li", "section"):
                flush()
                walk(c.children, fmt)
                flush()
                continue
            f = dict(fmt)
            if name in ("b", "strong"):
                f["bold"] = True
            elif name in ("i", "em"):
                f["italic"] = True
            elif name == "code":
                f["code"] = True
            elif name == "span":
                if "tag" in cls:
                    f["bold"] = True
                    f["color"] = next((TAG[k] for k in cls if k in TAG), None)
                col = css_color(style_of(c).get("color"))
                if col:
                    f["color"] = col
            walk(c.children, f)
            if name == "span" and "tag" in cls:
                cur.append(Run(" "))

    walk(nodes, {})
    flush()
    return out


def n_lines(runs: list, size: float, width_pt: float) -> int:
    words: list[list] = [[]]
    for r in runs:
        for i, part in enumerate(r.text.split(" ")):
            if i:
                words.append([])
            if part:
                words[-1].append((part, r))
    words = [w for w in words if w]
    space = text_w(" ", size)
    lines, x = 1, 0.0
    width_pt *= 0.93                       # margin for renderer differences
    for w in words:
        ww = sum(text_w(t, size, r.bold, r.italic, r.code) for t, r in w)
        if x == 0:
            x = ww
        elif x + space + ww <= width_pt:
            x += space + ww
        else:
            lines += 1
            x = ww
        while x > width_pt:
            lines += 1
            x -= width_pt
    return lines


BULLET_IN = 0.32


def para_h(p: Para, size: float, width_in: float) -> float:
    w = width_in - (BULLET_IN if p.kind == "li" else 0)
    runs = [Run(r.text, r.bold or p.kind == "h3", r.italic, r.code) for r in p.runs]
    return (n_lines(runs, size, w * 72) * size * LH + SPACE_AFTER) / 72


def paras_h(paras, size, width_in) -> float:
    return sum(para_h(p, size, width_in) for p in paras)


def single_line_w(p: Para, size: float) -> float:
    w = sum(text_w(r.text, size, r.bold or p.kind == "h3", r.italic, r.code) for r in p.runs)
    return w / 72 + (BULLET_IN if p.kind == "li" else 0)


def longest_word_w(p: Para, size: float) -> float:
    best = 0.0
    for r in p.runs:
        for t in r.text.split():
            best = max(best, text_w(t, size, r.bold or p.kind == "h3", r.italic, r.code))
    return best / 72


# --------------------------------------------------------------------------- pptx helpers

def pptx():
    from pptx import Presentation
    from pptx.dml.color import RGBColor
    from pptx.enum.shapes import MSO_SHAPE
    from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
    from pptx.oxml.ns import qn
    from pptx.util import Emu, Inches, Pt
    return dict(Presentation=Presentation, RGBColor=RGBColor, MSO_SHAPE=MSO_SHAPE,
                MSO_ANCHOR=MSO_ANCHOR, PP_ALIGN=PP_ALIGN, qn=qn, Emu=Emu, Inches=Inches, Pt=Pt)


P = None  # filled in main()


def set_bullet(para, ordered: bool):
    from lxml import etree
    qn, Inches = P["qn"], P["Inches"]
    pPr = para._p.get_or_add_pPr()
    pPr.set("marL", str(int(Inches(BULLET_IN))))
    pPr.set("indent", str(-int(Inches(BULLET_IN - 0.02))))
    clr = etree.SubElement(pPr, qn("a:buClr"))
    etree.SubElement(clr, qn("a:srgbClr")).set("val", BLUE)
    etree.SubElement(pPr, qn("a:buFont")).set("typeface", FONT)
    if ordered:
        etree.SubElement(pPr, qn("a:buAutoNum")).set("type", "arabicPeriod")
    else:
        etree.SubElement(pPr, qn("a:buChar")).set("char", "•")


def fill_tf(tf, paras, size, color=INK, h3_color=BLUE, align="l", insets=(0, 0, 0, 0)):
    Pt, Inches, RGB, AL, AN = P["Pt"], P["Inches"], P["RGBColor"], P["PP_ALIGN"], P["MSO_ANCHOR"]
    tf.word_wrap = True
    tf.auto_size = None
    tf.vertical_anchor = AN.TOP
    tf.margin_left, tf.margin_top, tf.margin_right, tf.margin_bottom = (Inches(v) for v in insets)
    for i, p in enumerate(paras):
        para = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        para.alignment = {"l": AL.LEFT, "c": AL.CENTER, "r": AL.RIGHT}[align]
        para.space_after = Pt(SPACE_AFTER)
        for r in p.runs:
            run = para.add_run()
            run.text = r.text
            f = run.font
            f.size = Pt(size)
            f.name = MONO if r.code else FONT
            f.bold = bool(r.bold or p.kind == "h3")
            f.italic = bool(r.italic)
            f.color.rgb = RGB.from_string(r.color or (h3_color if p.kind == "h3" else color))
        if p.kind == "li":
            set_bullet(para, p.ordered)


def add_text(slide, x, y, w, h, paras, size, **kw):
    Inches = P["Inches"]
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    fill_tf(tb.text_frame, paras, size, **kw)
    return tb


def add_rect(slide, x, y, w, h, fill, line=None, rounded=True, dash=False, line_w=1.0):
    Inches, Pt, RGB, SH = P["Inches"], P["Pt"], P["RGBColor"], P["MSO_SHAPE"]
    shp = slide.shapes.add_shape(SH.ROUNDED_RECTANGLE if rounded else SH.RECTANGLE,
                                 Inches(x), Inches(y), Inches(w), Inches(h))
    if rounded:
        shp.adjustments[0] = min(0.5, 0.08 / max(min(w, h), 0.01))
    if fill:
        shp.fill.solid()
        shp.fill.fore_color.rgb = RGB.from_string(fill)
    else:
        shp.fill.background()
    if line:
        shp.line.color.rgb = RGB.from_string(line)
        shp.line.width = Pt(line_w)
        if dash:
            from pptx.enum.dml import MSO_LINE
            shp.line.dash_style = MSO_LINE.DASH
    else:
        shp.line.fill.background()
    style = shp._element.find(P["qn"]("p:style"))
    if style is not None:                  # drop the theme style: no shadow, no theme colours
        shp._element.remove(style)
    return shp


# --------------------------------------------------------------------------- blocks

@dataclass
class Sizes:
    body: float = 18
    tdelta: int = 0
    cap: float = 99.0


class Block:
    def height(self, w, S) -> float: ...
    def place(self, slide, x, y, w, h, S) -> None: ...
    def images(self):
        return []


class TextBlock(Block):
    def __init__(self, paras, small=False, align="l"):
        self.paras, self.small, self.align = paras, small, align

    def size(self, S):
        return SMALL_PT if self.small else S.body

    def height(self, w, S):
        return paras_h(self.paras, self.size(S), w)

    def place(self, slide, x, y, w, h, S):
        add_text(slide, x, y, w, max(h, 0.2), self.paras, self.size(S),
                 color=MUTED if self.small else INK, align=self.align)


class BoxBlock(Block):
    PX, PY = 0.14, 0.09

    def __init__(self, paras, tone=None):
        self.paras, self.tone = paras, tone

    def height(self, w, S):
        return paras_h(self.paras, S.body, w - 2 * self.PX) + 2 * self.PY

    def place(self, slide, x, y, w, h, S):
        fill, line = BOX[self.tone]
        shp = add_rect(slide, x, y, w, h, fill, line)
        fill_tf(shp.text_frame, self.paras, S.body,
                insets=(self.PX, self.PY, self.PX, self.PY))


class ImageBlock(Block):
    def __init__(self, path: Path):
        from PIL import Image
        self.path = path
        with Image.open(path) as im:
            self.aspect = im.width / im.height

    def height(self, w, S):
        return min(w / self.aspect, S.cap)

    def place(self, slide, x, y, w, h, S):
        Inches = P["Inches"]
        iw = min(w, h * self.aspect)
        ih = iw / self.aspect
        slide.shapes.add_picture(str(self.path), Inches(x + (w - iw) / 2), Inches(y),
                                 Inches(iw), Inches(ih))

    def images(self):
        return [self]


class LegendBlock(Block):
    H = 0.34

    def __init__(self, items):
        self.items = items        # (text, fill, line, dashed)

    def height(self, w, S):
        return self.H

    def place(self, slide, x, y, w, h, S):
        widths = [text_w(t, SMALL_PT) / 72 + 0.3 for t, *_ in self.items]
        gap = 0.12
        cx = x + (w - sum(widths) - gap * (len(widths) - 1)) / 2
        for (t, fill, line, dashed), iw in zip(self.items, widths):
            shp = add_rect(slide, cx, y, iw, self.H, fill, line, dash=dashed)
            fill_tf(shp.text_frame, [Para([Run(t)])], SMALL_PT, align="c", insets=(0.05, 0.03, 0.05, 0.03))
            shp.text_frame.vertical_anchor = P["MSO_ANCHOR"].MIDDLE
            cx += iw + gap


class StripBlock(Block):
    H, SIZE = 0.4, 16

    def __init__(self, items):
        self.items = items        # list of Para

    def height(self, w, S):
        return self.H

    def place(self, slide, x, y, w, h, S):
        cx = x
        for p in self.items:
            for r in p.runs:
                if r.bold:
                    r.color = BLUE
            iw = single_line_w(p, self.SIZE) + 0.32
            shp = add_rect(slide, cx, y, iw, self.H, WHITE, RULE)
            fill_tf(shp.text_frame, [p], self.SIZE, align="c", insets=(0.08, 0.03, 0.08, 0.03))
            shp.text_frame.vertical_anchor = P["MSO_ANCHOR"].MIDDLE
            cx += iw + 0.12


class TableBlock(Block):
    PX, PY = 0.08, 0.045

    def __init__(self, el):
        cls = el.get("class") or []
        fs = style_of(el).get("font-size", "")
        self.base = 18 if "roomy" in cls else 16 if "mid" in cls else 14 if fs == "0.5em" else 15
        self.first_bold = "tools" in (el.find_parent("section").get("class") or []) or "tools" in cls
        self.rows = []            # list of [(paras, header, colspan)]
        for tr in el.find_all("tr"):
            row = []
            for td in tr.find_all(["th", "td"], recursive=False):
                paras = parse_nodes(td.children) or [Para([Run("")])]
                row.append((paras, td.name == "th", int(td.get("colspan", 1))))
            self.rows.append(row)
        self.ncols = max(sum(c[2] for c in r) for r in self.rows)

    def size(self, S):
        return max(12, self.base - S.tdelta)

    def _cells(self):
        for ri, row in enumerate(self.rows):
            ci = 0
            for paras, header, span in row:
                yield ri, ci, paras, header, span
                ci += span

    def _bold(self, paras, header, ci):
        if header or (self.first_bold and ci == 0):
            return [Para([Run(r.text, True, r.italic, r.code, r.color) for r in p.runs], p.kind) for p in paras]
        return paras

    def widths(self, w, size):
        pad = 2 * self.PX
        nat, minw = [0.4] * self.ncols, [0.4] * self.ncols
        for _, ci, paras, header, span in self._cells():
            if span != 1:
                continue
            paras = self._bold(paras, header, ci)
            nat[ci] = max(nat[ci], max(single_line_w(p, size) for p in paras) / 0.92 + pad)
            minw[ci] = max(minw[ci], max(longest_word_w(p, size) for p in paras) / 0.92 + pad)
        if self.first_bold:
            minw[0] = nat[0]
        self.too_narrow = sum(minw) > w       # a word would have to break: try a smaller font
        if sum(nat) <= w:
            return [v * w / sum(nat) for v in nat]
        lo, hi = 0.0, max(nat)
        for _ in range(60):
            t = (lo + hi) / 2
            if sum(max(m, min(n, t)) for n, m in zip(nat, minw)) > w:
                hi = t
            else:
                lo = t
        out = [max(m, min(n, lo)) for n, m in zip(nat, minw)]
        spare = w - sum(out)
        wide = [i for i, n in enumerate(nat) if n > lo]
        for i in wide:
            out[i] += spare / len(wide)
        return out

    def layout(self, w, S):
        size = self.size(S)
        cw = self.widths(w, size)
        rh = [0.0] * len(self.rows)
        for ri, ci, paras, header, span in self._cells():
            inner = sum(cw[ci:ci + span]) - 2 * self.PX
            rh[ri] = max(rh[ri], paras_h(self._bold(paras, header, ci), size, inner) + 2 * self.PY)
        return size, cw, rh

    def height(self, w, S):
        h = sum(self.layout(w, S)[2])
        return 99.0 if self.too_narrow and S.tdelta < 6 else h

    def place(self, slide, x, y, w, h, S):
        Inches, Emu, Pt, qn = P["Inches"], P["Emu"], P["Pt"], P["qn"]
        size, cw, rh = self.layout(w, S)
        frame = slide.shapes.add_table(len(self.rows), self.ncols, Inches(x), Inches(y),
                                       Inches(w), Inches(sum(rh)))
        tbl = frame.table
        tblPr = tbl._tbl.tblPr
        for attr in ("firstRow", "bandRow"):
            tblPr.set(attr, "0")
        for sid in tblPr.findall(qn("a:tableStyleId")):
            tblPr.remove(sid)
        emus = [int(Inches(v)) for v in cw]
        emus[-1] += int(Inches(w)) - sum(emus)
        for i, e in enumerate(emus):
            tbl.columns[i].width = Emu(e)
        for i, v in enumerate(rh):
            tbl.rows[i].height = Inches(v)
        for ri, ci, paras, header, span in self._cells():
            cell = tbl.cell(ri, ci)
            if span > 1:
                cell.merge(tbl.cell(ri, ci + span - 1))
            cell.margin_left = cell.margin_right = Inches(self.PX)
            cell.margin_top = cell.margin_bottom = Inches(self.PY)
            cell.vertical_anchor = P["MSO_ANCHOR"].TOP
            tf = cell.text_frame
            tf.word_wrap = True
            fill_tf(tf, self._bold(paras, header, ci), size,
                    insets=(self.PX, self.PY, self.PX, self.PY))
            header_row = all(c[1] for c in self.rows[ri])
            for k in range(ci, ci + span):
                c = tbl.cell(ri, k)
                set_borders(c, BLUE if header_row else RULE, 2.0 if header_row else 0.75)
                c.fill.solid()
                c.fill.fore_color.rgb = P["RGBColor"].from_string(BLUE_BG if header_row else WHITE)


def set_borders(cell, color, width_pt):
    from lxml import etree
    qn = P["qn"]
    tcPr = cell._tc.get_or_add_tcPr()
    for tag in ("a:lnL", "a:lnR", "a:lnT", "a:lnB"):
        for old in tcPr.findall(qn(tag)):
            tcPr.remove(old)
    for i, tag in enumerate(("a:lnL", "a:lnR", "a:lnT", "a:lnB")):
        ln = etree.Element(qn(tag))
        if tag == "a:lnB":
            ln.set("w", str(int(width_pt * 12700)))
            sf = etree.SubElement(ln, qn("a:solidFill"))
            etree.SubElement(sf, qn("a:srgbClr")).set("val", color)
        else:
            ln.set("w", "0")
            etree.SubElement(ln, qn("a:noFill"))
        tcPr.insert(i, ln)


class Stack(Block):
    def __init__(self, children, gap=0.12):
        self.children, self.gap = children, gap

    def height(self, w, S):
        hs = [c.height(w, S) for c in self.children]
        return sum(hs) + self.gap * max(0, len(hs) - 1)

    def place(self, slide, x, y, w, h, S):
        for c in self.children:
            ch = c.height(w, S)
            c.place(slide, x, y, w, ch, S)
            y += ch + self.gap

    def images(self):
        return [i for c in self.children for i in c.images()]


class Grid(Block):
    def __init__(self, children, ncols, fracs=None, gap=0.25):
        self.children, self.ncols, self.gap = children, ncols, gap
        self.fracs = fracs or [1.0] * ncols

    def col_ws(self, w):
        avail = w - self.gap * (self.ncols - 1)
        return [avail * f / sum(self.fracs) for f in self.fracs]

    def rows(self):
        return [self.children[i:i + self.ncols] for i in range(0, len(self.children), self.ncols)]

    def row_hs(self, w, S):
        cws = self.col_ws(w)
        return [max(c.height(cw, S) for c, cw in zip(row, cws)) for row in self.rows()]

    def height(self, w, S):
        hs = self.row_hs(w, S)
        return sum(hs) + self.gap * (len(hs) - 1)

    def place(self, slide, x, y, w, h, S):
        cws = self.col_ws(w)
        for row, rh in zip(self.rows(), self.row_hs(w, S)):
            cx = x
            for c, cw in zip(row, cws):
                c.place(slide, cx, y, cw, rh if isinstance(c, BoxBlock) else c.height(cw, S), S)
                cx += cw + self.gap
            y += rh + self.gap

    def images(self):
        return [i for c in self.children for i in c.images()]


GAPS = {"cols": 0.29, "cols3": 0.19, "cols4": 0.15, "card": 0.13}


def png_for(src: str) -> Path:
    path = PNG / (Path(src).stem + ".png")
    if not path.exists():
        sys.exit(f"missing {path}: run with --diagrams first")
    return path


def to_block(el) -> Block | None:
    name, cls = el.name, (el.get("class") or [])
    st = style_of(el)
    if name == "img":
        return ImageBlock(png_for(el["src"]))
    if name == "table":
        return TableBlock(el)
    if name in ("ul", "ol", "h3", "p"):
        paras = parse_nodes([el])
        if not paras:
            return None
        return TextBlock(paras, small="small" in cls,
                         align="c" if st.get("text-align") == "center" else "l")
    if name != "div":
        return None
    if "legend" in cls:
        items = []
        for sp in el.find_all("span", recursive=False):
            s = style_of(sp)
            items.append((sp.get_text(" ", strip=True), css_color(s.get("background")),
                          css_color(s.get("border-color")) or "9AA5B8",
                          "dashed" in s.get("border-style", "")))
        return LegendBlock(items)
    if "strip" in cls:
        return StripBlock([finish(parse_nodes([sp])[0].runs) for sp in el.find_all("span", recursive=False)])
    if "box" in cls:
        tone = next((t for t in ("blue", "green", "amber", "red") if t in cls), None)
        return BoxBlock(parse_nodes(el.children), tone)
    kind = next((k for k in ("cols4", "cols3", "card", "cols") if k in cls), None)
    kids = [b for b in (to_block(c) for c in el.find_all(recursive=False)) if b]
    if kind:
        ncols = {"cols": 2, "cols3": 3, "cols4": 4, "card": 2}[kind]
        if kind == "cols4":
            # Narrow boxes of the form <b>Label</b><br>text: put the label on the text's first line
            # (in blue) so the boxes are one line shorter and the diagram above them can be larger.
            for k in kids:
                if isinstance(k, BoxBlock) and len(k.paras) >= 2 and all(r.bold for r in k.paras[0].runs):
                    label = [Run(r.text, True, r.italic, r.code, BLUE) for r in k.paras[0].runs]
                    k.paras = [Para(label + [Run(" ")] + k.paras[1].runs)] + k.paras[2:]
        fracs = None
        m = st.get("grid-template-columns")
        if m:
            fracs = [float(v.replace("fr", "")) for v in m.split()]
            ncols = len(fracs)
        return Grid(kids, ncols, fracs, GAPS[kind])
    return Stack(kids) if kids else None


# --------------------------------------------------------------------------- slides

def notes_of(sec) -> str | None:
    a = sec.find("aside", class_="notes")
    return re.sub(r"\s+", " ", a.get_text(" ")).strip() if a else None


def footer(slide, sec, n, total, dark=False):
    src = sec.find("div", class_="src", recursive=False)
    color = "C9D1E0" if dark else MUTED
    if src:
        add_text(slide, X0, 7.05, 10.6, 0.3, [Para([Run(src.get_text(" ", strip=True))])], SRC_PT, color=color)
    add_text(slide, SLIDE_W - X0 - 1.2, 7.05, 1.2, 0.3, [Para([Run(f"{n} / {total}")])], SRC_PT,
             color=color, align="r")


def background(slide, color):
    fill = slide.background.fill
    fill.solid()
    fill.fore_color.rgb = P["RGBColor"].from_string(color)


def title_slide(slide, sec, n, total):
    background(slide, WHITE)
    add_rect(slide, 0, 0, 0.18, SLIDE_H, BLUE, rounded=False)
    kicker = sec.find(class_="kicker")
    y = 1.7
    add_text(slide, X0 + 0.2, y, W0 - 0.2, 0.4, parse_nodes([kicker]), 16, color=AMBER)
    h1 = parse_nodes([sec.find("h1")])[0]
    add_text(slide, X0 + 0.2, y + 0.5, W0 - 0.2, 1.0, [Para(h1.runs)], 48, color=INK)
    y += 1.75
    for p in sec.find_all("p", recursive=False):
        small = "small" in (p.get("class") or [])
        paras = parse_nodes([p])
        size = 16 if small else 22
        h = paras_h(paras, size, W0 - 0.2)
        add_text(slide, X0 + 0.2, y, W0 - 0.2, h, paras, size, color=MUTED if small else INK)
        y += h + 0.25
    footer(slide, sec, n, total)


def section_slide(slide, sec, n, total):
    bg = (css_color(sec.get("data-background-color")) or INK).upper()
    background(slide, bg)
    add_rect(slide, X0, 2.55, 1.1, 0.07, PEACH, rounded=False)
    kicker = sec.find(class_="kicker")
    add_text(slide, X0, 2.05, W0, 0.45, [Para([Run(kicker.get_text(strip=True).upper())])], 16,
             color=PEACH)
    h2 = parse_nodes([sec.find("h2")])[0]
    add_text(slide, X0, 2.8, W0, 1.0, [Para(h2.runs)], 40, color=WHITE)
    p = sec.find("p", recursive=False)
    if p:
        add_text(slide, X0, 3.85, W0 * 0.8, 1.2, parse_nodes([p]), 22, color="E3E8F1")
    footer(slide, sec, n, total, dark=True)


def content_slide(slide, sec, n, total):
    background(slide, WHITE)
    y = 0.32
    kicker = sec.find("div", class_="kicker", recursive=False)
    if kicker:
        add_text(slide, X0, y, W0, 0.3, [Para([Run(kicker.get_text(strip=True).upper(), bold=True)])],
                 KICKER_PT, color=AMBER)
        y += 0.34
    title = Para(parse_nodes([sec.find("h2", recursive=False)])[0].runs)
    th = n_lines([Run(r.text, True) for r in title.runs], TITLE_PT, W0 * 72) * TITLE_PT * 1.22 / 72
    add_text(slide, X0, y, W0, th, [title], TITLE_PT, color=INK, h3_color=INK)
    for r in slide.shapes[-1].text_frame.paragraphs[0].runs:
        r.font.bold = True
    y += th + 0.08
    add_rect(slide, X0, y, 0.9, 0.05, BLUE, rounded=False)
    top = y + 0.2
    skip = {"h2", "aside"}
    kids = [c for c in sec.find_all(recursive=False)
            if c.name not in skip and not ({"src", "kicker"} & set(c.get("class") or []))]
    body = Stack([b for b in (to_block(c) for c in kids) if b], gap=0.16)
    avail = 6.95 - top
    imgs = body.images()
    chosen = None
    steps = [(18, 0), (18, 1), (17, 1), (17, 2), (16, 2), (16, 3), (16, 4), (16, 5), (16, 6)]
    for b, t in steps:
        S = Sizes(body=b, tdelta=t, cap=MIN_IMG if imgs else 99)
        if body.height(W0, S) <= avail:
            chosen = S
            break
    if chosen is None:
        chosen = Sizes(body=16, tdelta=6, cap=MIN_IMG if imgs else 99)
        if not imgs:
            print(f"WARNING slide {n}: content {body.height(W0, chosen):.2f} in > {avail:.2f} in",
                  file=sys.stderr)
    if imgs:
        lo, hi = (MIN_IMG if body.height(W0, chosen) <= avail else 1.0), 6.0
        for _ in range(40):
            mid = (lo + hi) / 2
            chosen.cap = mid
            if body.height(W0, chosen) <= avail:
                lo = mid
            else:
                hi = mid
        chosen.cap = lo
    body.place(slide, X0, top, W0, avail, chosen)
    footer(slide, sec, n, total)
    return chosen


# --------------------------------------------------------------------------- update section

UPDATE_HTML = """
<section class="section-title" data-background-color="#14213d">
  <div class="kicker">Update</div>
  <h2>Update, 12 November 2026</h2>
  <p>Since this deck (2026-10-28): what is built, what a search gives, the decisions taken, and what comes next.</p>
  <div class="src">planning/tasks/STATUS.md (handoff 2026-11-12 10:30 CST; merged cards) · planning/DECISIONS.md DEC-24 to DEC-29</div>
  <aside class="notes">These five slides come only from planning/tasks/STATUS.md (the handoff of 2026-11-12 10:30 CST and the table of merged cards) and from planning/DECISIONS.md rows DEC-24 to DEC-29. Every number carries the label those files give it; where they give none, the number is a count. The slides before this section are the 2026-10-28 deck, unchanged.</aside>
</section>

<section>
  <h2>95 cards are merged and real releases verify on the pod; the 50 Gbp release rM.1 has its clusters but stopped at the token step.</h2>
  <div class="cols">
    <div class="box green"><h3>Built and verified</h3><ul>
      <li><code>plan/task-cards</code>: 95 cards merged; last green <code>ci</code> and <code>cross-host</code> at <code>3ff6bc4</code> (2026-11-12)</li>
      <li><code>main</code> has pull requests 1 and 2 (87 cards) and is behind</li>
      <li>First real <code>build-stage</code> run (C-10): 1,947 phages and 250 prokaryotes, two releases, both verify</li>
      <li>Releases <code>rN.1a</code>, <code>rN.1b</code>, <code>rN.1</code>, <code>rN.2</code>, <code>rN.3</code> all verify (C-19)</li>
      <li>Search pack of <code>rN.1</code>: 5.51 GB, built at 750 MiB peak (<i>measured</i>, W-16)</li>
      <li>Data root holds stages S1 and P1 to P6, on local disk, not yet copied to <code>/stor</code></li>
    </ul></div>
    <div class="box amber"><h3>The 50 Gbp build (unit <code>r50-driver</code>)</h3><ul>
      <li>Stages P3 to P6 ingested and gene-called</li>
      <li><code>rM.1</code> published its clustering at 10:53 CST: 45,815,843 proteins, 31,183,876 NR90 and 8,192,065 NR50 clusters (<i>measured</i>)</li>
      <li><code>rM.1.tokens</code> failed at 10:57: DuckDB out of memory in the slot join (card B-46)</li>
      <li>Relaunched with <code>SKIP_RM1=1</code>: Pfam scan of P3 to P6 to about midday Saturday (<i>estimate</i>), then <code>rM.2</code> once B-46 is merged</li>
    </ul></div>
  </div>
  <div class="src">planning/tasks/STATUS.md: handoff 2026-11-12 10:30 CST; merged cards C-10, C-19, W-16</div>
  <aside class="notes">What exists on 12 November. The 10 Gbp-era releases rN.1a to rN.3 all verify; the 50 Gbp release rM.1 has clusters but no tokens yet, because the token step ran out of memory (card B-46). The driver carries on with the Pfam scan meanwhile. Note the data root is on local disk that a reboot clears.</aside>
</section>

<section>
  <h2>A search runs end to end on a real release: one command gives hits, sequences and neighbour summaries, and the worker peak fell from 26.27 GiB to 1.78 GiB.</h2>
  <div class="cols">
    <div>
      <h3>What a user runs and gets</h3>
      <ul>
        <li><code>bio-box search hmm --hmm FILE --release R --bits 30 --neighbors 5 --positive ACCESSION:START-END --out DIR</code> (W-07)</li>
        <li>Recipe, hit table, FASTA, neighbour summary, per-locus tables, and a <code>README.txt</code> naming every column (W-07, W-20)</li>
        <li>0-based half-open coordinates; strand <code>+1</code>/<code>-1</code> in every exported file (W-20)</li>
        <li><code>bio-box releases</code>, <code>seq --protein UID</code>, <code>seq --contig ACC:START-END</code>, read-only <code>bio-box sql FILE</code> (W-23)</li>
      </ul>
    </div>
    <div>
      <h3>The trial search on <code>rN.1</code> (<i>measured</i>, wilkcomp02)</h3>
      <table class="mid">
        <tr><th></th><th>Before W-19</th><th>After W-19</th><th>After W-21</th></tr>
        <tr><td>Worker peak</td><td>26.27 GiB</td><td>1.34 GiB</td><td>1.78 GiB</td></tr>
        <tr><td>Wall time</td><td>534 s</td><td>294 s</td><td>244.6 s</td></tr>
        <tr><td>Proteins found</td><td>1,957</td><td>1,957</td><td>1,957</td></tr>
      </table>
      <p class="small">W-21 ran on a loaded pod; its peak rose 33% over W-19, cause not fully established. Also <i>measured</i> on <code>rN.1</code>: <code>partners_of</code> with 3,509 anchors 12.35 s before, 0.52 s after (B-38). Not yet: no tiering step exists (W-20).</p>
    </div>
  </div>
  <div class="src">planning/tasks/STATUS.md merged cards W-07, W-19, W-20, W-21, W-23, B-38</div>
  <aside class="notes">The search path that the deck promised now runs on a real release. W-19 made the search stream its target, which took the worker from 26.27 GiB to 1.34 GiB with byte-identical output; W-21 added read-ahead, faster but 33% more memory. Every number here is measured on wilkcomp02 and comes from the merged-cards table.</aside>
</section>

<section>
  <h2>Six decisions since 10 November: each release clusters from scratch under its own CPU cap, and a contig may belong to several genomes.</h2>
  <table class="mid">
    <tr><th>Row</th><th>Date</th><th>Decision</th><th>Basis, as labelled in DECISIONS.md</th></tr>
    <tr><td>DEC-24</td><td>2026-11-10</td><td>Each release clusters all its proteins from scratch; <code>mmseqs clusterupdate</code> is not used</td><td>Our side of clustering peaks at about 1 GB at 2.05 M proteins (<i>measured</i>); about 36 h per release at 95 M proteins (<i>derived</i>)</td></tr>
    <tr><td>DEC-25</td><td>2026-11-10</td><td>Clustering is exempt from the 300 CPU-minute rule until it is split, under its own cap</td><td>10,620,684 proteins used 2,195.9 CPU-minutes in one process (<i>measured</i>)</td></tr>
    <tr><td>DEC-26</td><td>2026-11-11</td><td>A contig shared by stages belongs to the stage that stored it first</td><td>4 of 54,265 contigs of the 50 Gbp set are shared across stages (<i>measured</i>)</td></tr>
    <tr><td>DEC-27</td><td>2026-11-11</td><td>A contig may belong to several genomes (amends DEC-26); repeated records are no longer dropped</td><td>The input script had removed 6,818 of 34,076 INPHARED records from stage S1 (<i>measured</i>)</td></tr>
    <tr><td>DEC-28</td><td>2026-11-11</td><td>ncRNA: 588 curated Rfam 15.1 models on intergenic regions, with clan competition</td><td>Full Rfam 1,910 CPU-hours per 10 Gbp; this set 247 at 99.4% of hits (<i>derived</i> from per-model CPU, <i>measured</i>)</td></tr>
    <tr><td>DEC-29</td><td>2026-11-12</td><td>A clustering call that cannot be split under the budget runs as one process, under clustering's cap</td><td>Starts at about 74 M proteins (NR90 step) and 108 M (NR50 step 0) (<i>derived</i>)</td></tr>
  </table>
  <p class="small">DEC-24 replaces the incremental <code>clusterupdate</code> path shown earlier in this deck. All six were decided by the PI.</p>
  <div class="src">planning/DECISIONS.md rows DEC-24 to DEC-29</div>
  <aside class="notes">Two of these change what the deck said: clusters are no longer updated incrementally (DEC-24), and genome membership is many-to-many (DEC-27). DEC-25 and DEC-29 are the CPU-budget exemptions that from-scratch clustering needs on the pod. DEC-28 cut the ncRNA scan from full Rfam to a curated set.</aside>
</section>

<section>
  <h2>Next: fix the token step (B-46) so rM.2 can be built, copy the data root to permanent storage, then add huge phages and PLSDB.</h2>
  <div class="cols" style="grid-template-columns: 1.4fr 1fr">
    <div>
      <h3>Next, in order</h3>
      <ol>
        <li>B-46 first: it blocks <code>rM.2</code>. Merge C-21b, W-22, C-23, B-45 as they report; then C-21c, C-16, W-24, W-15</li>
        <li>When <code>rM.2</code> is built and verified: copy the data root to <code>/stor</code> (one <code>rsync</code>, atomic rename); clear obsolete scratch</li>
        <li>C-24 (stage H1, huge phages) after the Pfam scan of P3 to P6; B-42 after B-45; then PLSDB (no card yet; needs C-21b and C-21c)</li>
        <li>A pull request bringing <code>main</code> up to <code>plan/task-cards</code>; runbook and <code>pod_pytest.sh</code> on the new image</li>
        <li>Finish measurement 0-08a (Rfam sample scan); S2C phase 2 (card C-25) waits for the PI</li>
      </ol>
    </div>
    <div>
      <div class="box blue"><h3>Agents running</h3>C-21b genome membership view; W-22 one pack directory per host; C-23 run records name their code; B-45 unsplittable clustering call as one process (DEC-29).</div>
      <div class="box amber" style="margin-top:10px"><h3>On the pod</h3>Pfam scan of P3 to P6, to about midday Saturday (<i>estimate</i>). Rfam sample scan 0-08a: 16 shards done at 10:12 CST. The data root is on local disk, which a reboot clears.</div>
    </div>
  </div>
  <div class="src">planning/tasks/STATUS.md: handoff 2026-11-12 10:30 CST, "Agents running", "On the pod now", "Next, in order"</div>
  <aside class="notes">The order matters: B-46 unblocks rM.2, and only a verified rM.2 justifies the copy to /stor. Huge phages (H1) and PLSDB come after; PLSDB needs the genome-membership cards C-21b and C-21c first. Items waiting for the PI are listed in STATUS.md under "Open items for the PI".</aside>
</section>
"""


# --------------------------------------------------------------------------- main

def build() -> None:
    global P
    P = pptx()
    from bs4 import BeautifulSoup

    soup = BeautifulSoup(HTML.read_text(encoding="utf-8"), "html.parser")
    sections = soup.select("div.slides > section")
    update = BeautifulSoup(UPDATE_HTML, "html.parser").find_all("section", recursive=False)
    appendix = next(i for i, s in enumerate(sections)
                    if "section-title" in (s.get("class") or [])
                    and s.find(class_="kicker").get_text(strip=True) == "Appendix")
    deck = sections[:appendix] + update + sections[appendix:]

    prs = P["Presentation"]()
    prs.slide_width, prs.slide_height = P["Inches"](SLIDE_W), P["Inches"](SLIDE_H)
    layout = prs.slide_layouts[6]
    total = len(deck)
    for n, sec in enumerate(deck, 1):
        slide = prs.slides.add_slide(layout)
        if sec.find("h1"):
            title_slide(slide, sec, n, total)
        elif "section-title" in (sec.get("class") or []):
            section_slide(slide, sec, n, total)
        else:
            S = content_slide(slide, sec, n, total)
            print(f"slide {n:2d}: body {S.body} pt, table delta {S.tdelta}"
                  + (f", diagram cap {S.cap:.2f} in" if S.cap < 99 else ""))
        notes = notes_of(sec)
        if notes:
            slide.notes_slide.notes_text_frame.text = notes

    cp = prs.core_properties
    cp.title = "bio-box staged build"
    cp.author = cp.last_modified_by = "bio-box"
    cp.created = cp.modified = dt.datetime(2026, 11, 12)
    cp.revision = 1
    prs.save(OUT)
    print(f"wrote {OUT.relative_to(HERE.parent)}: {total} slides "
          f"({len(sections)} from the HTML deck, {len(update)} update)")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--diagrams", action="store_true",
                    help="render diagrams/d*.svg to diagrams/png/ with Chrome, then stop")
    args = ap.parse_args()
    if args.diagrams:
        render_diagrams()
    else:
        build()


if __name__ == "__main__":
    main()
