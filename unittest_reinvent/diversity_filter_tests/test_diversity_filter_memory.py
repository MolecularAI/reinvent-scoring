import pandas as pd

from unittest_reinvent.diversity_filter_tests.test_diversity_filter_memory_base import BaseDiversityFilterMemory


class TestDiversityFilterMemoryAdd(BaseDiversityFilterMemory):

    def test_empty(self):
        self.assertFalse(self.mem.smiles_exists(self.smi))
        self.assertEqual(self.mem.number_of_smiles(), 0)

    def test_one_add(self):
        self.mem._add_to_memory_dataframe(0, self.smi, self.scaffold, {})

        self.assertTrue(self.mem.smiles_exists(self.smi))
        self.assertFalse(self.mem.smiles_exists(self.smi2))
        self.assertEqual(self.mem.number_of_smiles(), 1)
        self.assertEqual(self.mem.number_of_scaffolds(), 1)
        self.assertEqual(self.mem.scaffold_instances_count(self.scaffold), 1)
        self.assertEqual(len(self.mem.get_memory()), 1)

    def test_two_adds(self):
        self.mem._add_to_memory_dataframe(0, self.smi, self.scaffold, {})
        self.mem._add_to_memory_dataframe(0, self.smi2, self.scaffold, {})

        self.assertTrue(self.mem.smiles_exists(self.smi))
        self.assertTrue(self.mem.smiles_exists(self.smi2))
        self.assertEqual(self.mem.number_of_smiles(), 2)
        self.assertEqual(self.mem.number_of_scaffolds(), 1)
        self.assertEqual(self.mem.scaffold_instances_count(self.scaffold), 2)
        self.assertEqual(len(self.mem.get_memory()), 2)

    def test_custom_component_scores(self):
        self.mem._add_to_memory_dataframe(self.step, self.smi, self.scaffold, {"foo": 2, "bar": 3})
        df = self.mem.get_memory()

        self.assertEqual(df.loc[0, "foo"], 2)
        self.assertEqual(df.loc[0, "bar"], 3)
        self.assertTrue(
            all(c in df.columns for c in {"foo", "bar", "SMILES", "Scaffold", "Step"})
        )

    def test_set_memory(self):
        df = pd.DataFrame({"x": [1, 2, 3]})
        self.mem.set_memory(df)

        self.assertEqual(len(self.mem.get_memory()), len(df))

    def test_add_duplicate_smi(self):
        self.mem._add_to_memory_dataframe(self.step, self.smi, self.scaffold, {})
        self.mem._add_to_memory_dataframe(self.step+1, self.smi, self.scaffold, {})

        self.assertTrue(self.mem.smiles_exists(self.smi))
        self.assertEqual(self.mem.number_of_smiles(), 1)
        self.assertEqual(self.mem.get_memory().loc[0, "Step"], self.step, "Expected original step.")
