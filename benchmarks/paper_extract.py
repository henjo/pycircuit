# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Pull figures and tables out of the reference papers as images.

Why this exists: the tables in the DDD papers are **300-dpi bitonal scans**, not
text.  ``pdftotext`` returns their captions and nothing else, so the reference
numbers we want to calibrate against -- matrix sizes, vertex counts, the
partitions used -- are unreachable by any text tool.  Rendering the relevant page
or embedded image to PNG makes them readable.

This matters beyond one table.  Several papers in ``doc/ddd_references.md`` carry
benchmark circuits and published figures we may want to check our own
measurements against, and all of them are scans of the same vintage.

Usage::

    python benchmarks/paper_extract.py list   PAPER.pdf
    python benchmarks/paper_extract.py page   PAPER.pdf 11
    python benchmarks/paper_extract.py images PAPER.pdf 11
    python benchmarks/paper_extract.py find   PAPER.pdf "TABLE IV"

``find`` reports which page a caption appears on, which is usually the quickest
way to locate a table worth rendering.

Requires ``pdftoppm``, ``pdfimages`` and ``pdftotext`` (poppler-utils).
"""

import argparse
import os
import re
import subprocess
import sys


def _run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=False)


def _need(tool):
    if _run(['which', tool]).returncode != 0:
        sys.exit('%s not found; install poppler-utils' % tool)


def page_count(pdf):
    out = _run(['pdfinfo', pdf]).stdout
    match = re.search(r'^Pages:\s+(\d+)', out, re.M)
    return int(match.group(1)) if match else 0


def cmd_list(args):
    """Show page count and where the embedded images are."""
    _need('pdfimages')
    print('%s -- %d pages' % (os.path.basename(args.pdf), page_count(args.pdf)))
    out = _run(['pdfimages', '-list', args.pdf]).stdout
    lines = out.splitlines()
    print('\n'.join(lines[:2]))
    for line in lines[2:]:
        print(line)


def cmd_find(args):
    """Report which pages contain a piece of text -- usually a caption."""
    _need('pdftotext')
    needle = args.text.lower()
    hits = []
    for page in range(1, page_count(args.pdf) + 1):
        text = _run(['pdftotext', '-f', str(page), '-l', str(page),
                     args.pdf, '-']).stdout
        if needle in text.lower():
            hits.append(page)
    if not hits:
        print('not found: %r' % args.text)
        return
    print('%r appears on page(s): %s' % (args.text,
                                         ', '.join(map(str, hits))))
    print('render one with:  page %s <n>' % args.pdf)


def cmd_page(args):
    """Render whole pages to PNG.

    Rendering the page rather than extracting the embedded image keeps the
    caption and column headings attached, which is what makes a table
    interpretable.
    """
    _need('pdftoppm')
    os.makedirs(args.outdir, exist_ok=True)
    stem = os.path.join(args.outdir,
                        os.path.splitext(os.path.basename(args.pdf))[0])
    for page in args.pages:
        prefix = '%s_p%d' % (stem, page)
        result = _run(['pdftoppm', '-png', '-r', str(args.dpi),
                       '-f', str(page), '-l', str(page), args.pdf, prefix])
        if result.returncode != 0:
            print('page %d failed: %s' % (page, result.stderr.strip()))
            continue
        for name in sorted(os.listdir(args.outdir)):
            if name.startswith(os.path.basename(prefix)):
                print(os.path.join(args.outdir, name))


def cmd_images(args):
    """Extract the embedded images of a page -- tighter crops than a full page."""
    _need('pdfimages')
    os.makedirs(args.outdir, exist_ok=True)
    stem = os.path.join(args.outdir,
                        os.path.splitext(os.path.basename(args.pdf))[0])
    for page in args.pages:
        prefix = '%s_p%d_img' % (stem, page)
        result = _run(['pdfimages', '-png', '-f', str(page), '-l', str(page),
                       args.pdf, prefix])
        if result.returncode != 0:
            print('page %d failed: %s' % (page, result.stderr.strip()))
            continue
        for name in sorted(os.listdir(args.outdir)):
            if name.startswith(os.path.basename(prefix)):
                path = os.path.join(args.outdir, name)
                print('%s  (%d bytes)' % (path, os.path.getsize(path)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest='command', required=True)

    p = sub.add_parser('list', help='page count and embedded images')
    p.add_argument('pdf')
    p.set_defaults(func=cmd_list)

    p = sub.add_parser('find', help='locate a caption')
    p.add_argument('pdf')
    p.add_argument('text')
    p.set_defaults(func=cmd_find)

    p = sub.add_parser('page', help='render pages to PNG')
    p.add_argument('pdf')
    p.add_argument('pages', nargs='+', type=int)
    p.add_argument('-r', '--dpi', type=int, default=200)
    p.add_argument('-o', '--outdir', default='paper_figures')
    p.set_defaults(func=cmd_page)

    p = sub.add_parser('images', help='extract a page\'s embedded images')
    p.add_argument('pdf')
    p.add_argument('pages', nargs='+', type=int)
    p.add_argument('-o', '--outdir', default='paper_figures')
    p.set_defaults(func=cmd_images)

    args = ap.parse_args(argv)
    return args.func(args) or 0


if __name__ == '__main__':
    sys.exit(main())
