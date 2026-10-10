"""The workbench UI must stay bare: no suggested structures or chemicals.

A product requirement, not a style preference. What a researcher chooses to
analyze should not be nudged by the tool — a placeholder reading "aspirin,
ibuprofen, caffeine" or a one-click "1HSG" chip is a suggestion, and
suggestions belong on the help page where they read as documentation rather
than as a prompt.

These tests guard the workbench surface (``index.html`` and ``app.js``).
``help.html`` is deliberately exempt: documenting reference cases is its job.
"""

import os
import re
import unittest

WEB = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "web")

# Workbench surface — what a user sees before they have chosen anything.
WORKBENCH = ("index.html", "app.js")

# Names that only appear in this codebase as example molecules. Deliberately a
# small, concrete list rather than a clever heuristic: it is readable, it does
# not produce false positives on scientific prose, and adding to it is obvious.
EXAMPLE_MOLECULES = (
    "aspirin", "ibuprofen", "caffeine", "benzamidine", "indinavir",
    "methotrexate", "paracetamol", "acetaminophen", "penicillin", "warfarin",
)

# PDB entries that have historically been used here as demo structures.
EXAMPLE_PDB_IDS = ("1HSG", "1CA2", "4HHB", "3PTB", "4DFR", "6BCX")


def _read(name):
    with open(os.path.join(WEB, name), encoding="utf-8") as fh:
        return fh.read()


class TestWorkbenchIsBare(unittest.TestCase):
    def test_no_example_molecule_names(self):
        for name in WORKBENCH:
            text = _read(name).lower()
            for mol in EXAMPLE_MOLECULES:
                self.assertNotIn(
                    mol, text,
                    f"{name} names the example molecule '{mol}'. The workbench "
                    f"must not suggest what to analyze; put it in help.html.",
                )

    def test_no_example_pdb_ids(self):
        for name in WORKBENCH:
            text = _read(name)
            for pdb in EXAMPLE_PDB_IDS:
                self.assertNotIn(
                    pdb, text,
                    f"{name} names the example structure '{pdb}'. "
                    f"Put reference entries in help.html.",
                )

    def test_placeholders_describe_the_field_not_a_value(self):
        """A placeholder may say what kind of value a field takes."""
        html = _read("index.html")
        placeholders = re.findall(r'placeholder="([^"]*)"', html)
        self.assertTrue(placeholders, "expected some placeholders to check")
        for ph in placeholders:
            low = ph.lower()
            for mol in EXAMPLE_MOLECULES:
                self.assertNotIn(mol, low, f"placeholder suggests a molecule: {ph!r}")
            # A bare number in a placeholder is almost always a sample CID.
            self.assertNotRegex(
                ph, r"\b\d{3,}\b", f"placeholder looks like a sample id: {ph!r}"
            )

    def test_no_prefilled_input_values(self):
        """An input carrying value= would hand the user a starting molecule."""
        html = _read("index.html")
        for tag in re.findall(r"<input\b[^>]*>", html, re.I):
            if 'type="checkbox"' in tag.lower() or 'type="file"' in tag.lower():
                continue
            self.assertNotIn(
                "value=", tag.lower(), f"input ships with a prefilled value: {tag}"
            )

    def test_no_try_this_phrasing(self):
        for name in WORKBENCH:
            text = _read(name)
            for pattern in (r"\bTry\s+\d", r"\be\.g\.\s*\d", r"\bfor example,?\s*\d"):
                self.assertNotRegex(
                    text, pattern,
                    f"{name} suggests a specific entry to try; move it to help.html",
                )


class TestHelpPageCarriesTheGuidance(unittest.TestCase):
    """Whatever the workbench stopped offering has to be documented somewhere."""

    def test_help_page_exists_and_is_linked_everywhere(self):
        self.assertTrue(os.path.isfile(os.path.join(WEB, "help.html")))
        for name in ("index.html", "privacy.html", "terms.html", "api.html"):
            self.assertIn("/help.html", _read(name), f"{name} does not link to help")

    def test_help_covers_the_topics_the_workbench_no_longer_explains(self):
        html = _read("help.html").lower()
        for topic in ("workflow", "glossary", "limitation", "data source",
                      "benchmark", "panel"):
            self.assertIn(topic, html, f"help page does not cover '{topic}'")

    def test_help_does_not_hardcode_the_reference_cases(self):
        """They come from /api/benchmark/cases, so the two cannot drift apart."""
        html = _read("help.html")
        for pdb in EXAMPLE_PDB_IDS:
            self.assertNotIn(pdb, html)
        self.assertIn("benchCases", html)

    def test_help_states_that_nothing_is_preloaded(self):
        # Whitespace-normalised: the sentence is line-wrapped in the source.
        text = " ".join(_read("help.html").lower().split())
        self.assertIn("no suggested chemicals", text)
        self.assertIn("nothing is pre-loaded", text)


if __name__ == "__main__":
    unittest.main()
