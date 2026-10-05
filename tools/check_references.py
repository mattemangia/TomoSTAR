#!/usr/bin/env python3
# Copyright 2026 Matteo Mangiagalli
# SPDX-License-Identifier: Apache-2.0
"""
Checks the references of a Markdown file against the registries of their DOIs.

Every item of the "## References" section that carries a DOI (https://doi.org/...) is resolved on
Crossref (journal articles, books, chapters) or, when Crossref does not know it, on DataCite (data
sets, software, FDSN network DOIs). The registered metadata are then compared with the reference
as written: the year, the surname of the first author (or the name of the first corporate creator)
and, when the registry has them, the volume, the issue and the first page. Items without a DOI are
listed so that they can be checked by hand (older articles and book chapters that were never
registered).

    python3 tools/check_references.py [README.md] [--mailto ADDRESS]

The exit status is 1 when a DOI does not resolve or its metadata disagree with the reference.
Only the Python standard library is needed.
"""
import argparse
import json
import re
import sys
import unicodedata
import urllib.error
import urllib.parse
import urllib.request


def fold(text):
    """Lower case without accents and with ASCII hyphens, for comparing names as registries and authors write them."""
    text = text.replace("\u2010", "-").replace("\u2011", "-")
    return "".join(c for c in unicodedata.normalize("NFKD", text) if not unicodedata.combining(c)).lower()


def get_json(url, mailto):
    request = urllib.request.Request(url, headers={"User-Agent": f"TomoSTAR-reference-check (mailto:{mailto})"})
    try:
        with urllib.request.urlopen(request, timeout=60) as response:
            return json.load(response)
    except urllib.error.HTTPError as error:
        if error.code == 404:
            return None
        raise


def crossref(doi, mailto):
    data = get_json("https://api.crossref.org/works/" + urllib.parse.quote(doi), mailto)
    if data is None:
        return None
    m = data["message"]
    # The year of the issue counts too: an article online in one year is cited with its issue's year.
    issue_print = (m.get("journal-issue") or {}).get("published-print") or {}
    years = {p[0][0] for p in [(m.get(k) or {}).get("date-parts") or [[None]] for k in ("published-print", "published-online", "issued")]
             + [issue_print.get("date-parts") or [[None]]] if p[0][0]}
    if not years:
        return None  # an incomplete record (a DOI of another agency): the caller asks DataCite
    authors = [a.get("family") or a.get("name", "") for a in m.get("author", [])] or \
              [a.get("family") or a.get("name", "") for a in m.get("editor", [])]
    return {"registry": "Crossref", "title": (m.get("title") or [""])[0], "years": years,
            "authors": authors, "volume": m.get("volume"), "issue": m.get("issue"),
            "first_page": (m.get("page") or "").split("-")[0] or None}


def datacite(doi, mailto):
    data = get_json("https://api.datacite.org/dois/" + urllib.parse.quote(doi), mailto)
    if data is None:
        return None
    a = data["data"]["attributes"]
    authors = [c.get("familyName") or c.get("name", "") for c in a.get("creators", [])]
    return {"registry": "DataCite", "title": a["titles"][0]["title"], "years": {int(a["publicationYear"])},
            "authors": authors, "volume": None, "issue": None, "first_page": None}


def references(markdown):
    section = markdown.split("## References", 1)[1]
    section = re.split(r"\n## ", section, maxsplit=1)[0]
    # Items of the list that carry a year in parentheses (the references, not other text of the section).
    return [re.sub(r"\s+", " ", item).strip() for item in re.split(r"\n- ", "\n" + section) if re.search(r"\(\d{4}\)", item)]


def check(item, mailto):
    """Returns (doi, registry, problems) for one reference; doi is None when it has none."""
    match = re.search(r"https://doi\.org/(\S+)", item)
    if not match:
        return None, None, []
    doi = match.group(1).rstrip(".")
    meta = crossref(doi, mailto) or datacite(doi, mailto)
    if meta is None:
        return doi, None, ["the DOI does not resolve on Crossref or DataCite"]
    problems = []
    year = int(re.search(r"\((\d{4})\)", item).group(1))
    if year not in meta["years"]:
        problems.append(f"year {year}, registered {sorted(meta['years'])}")
    first = fold(meta["authors"][0]) if meta["authors"] else ""
    lead = fold(item.split("(")[0])
    # The first registered creator must appear among the authors written before the year. Corporate
    # creators are compared on their first word (for example "Istituto" for INGV is accepted when the
    # reference writes the acronym, which the registry gives in parentheses).
    if first and first not in lead and first.split()[0] not in lead and not any(
            w.strip("()") in lead for w in first.split() if w.isupper() or w.startswith("(")):
        acronym = re.findall(r"\(([a-z\-]+)\)", first)
        if not any(a in lead for a in acronym):
            problems.append(f"first author '{meta['authors'][0]}' not in '{item.split('(')[0].strip()}'")
    after = item.split(")", 1)[1]
    vol = re.search(r", (\d+)(?:\((\w+)\))?, (\w+)(?:-\w+)?\.", after)
    if vol and meta["volume"]:
        volume, issue, page = vol.groups()
        if volume != str(meta["volume"]):
            problems.append(f"volume {volume}, registered {meta['volume']}")
        if issue and meta["issue"] and issue != str(meta["issue"]):
            problems.append(f"issue {issue}, registered {meta['issue']}")
        if meta["first_page"] and page.lower() != meta["first_page"].lower():
            problems.append(f"first page {page}, registered {meta['first_page']}")
    return doi, meta["registry"], problems


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("markdown", nargs="?", default="README.md")
    parser.add_argument("--mailto", default="tomostar@users.noreply.github.com",
                        help="contact address sent to the registries (their polite-pool etiquette)")
    args = parser.parse_args()
    failed = 0
    without = []
    items = references(open(args.markdown, encoding="utf-8").read())
    for item in items:
        doi, registry, problems = check(item, args.mailto)
        label = item[:60] + ("..." if len(item) > 60 else "")
        if doi is None:
            without.append(label)
            continue
        status = "OK  " if not problems else "FAIL"
        failed += bool(problems)
        print(f"{status} {registry or '-':8s} {doi:40s} {label}")
        for p in problems:
            print(f"       {p}")
    for label in without:
        print(f"--   no DOI   {'':40s} {label}")
    print(f"{len(items)} references, {len(items) - len(without)} with a DOI, {failed} with problems.")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
