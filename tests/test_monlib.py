"""
Author: "Keitaro Yamashita, Garib N. Murshudov"
MRC Laboratory of Molecular Biology
    
This software is released under the
Mozilla Public License, version 2.0; see LICENSE.
"""
from __future__ import absolute_import, division, print_function, generators
import unittest
import os
import glob
import gemmi

monlib_path = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

def read_mon_lib_list():
    doc = gemmi.cif.read(os.path.join(monlib_path, "list", "mon_lib_list.cif"))
    comp_list = doc.find_block("comp_list")
    link_list = doc.find_block("link_list")
    mod_list = doc.find_block("mod_list")
    return comp_list, link_list, mod_list


class TestMonlib(unittest.TestCase):
    def setUp(self): self.errors = []
    def tearDown(self): self.assertEqual(len(self.errors), 0, msg="\n"+"\n".join(self.errors))

    def test_list_consistency(self):
        comp_list, link_list, mod_list = read_mon_lib_list()
        lgroups = {row.str(0):row.str(1) for row in comp_list.find("_chem_comp.", ["id", "group"])}
        cgroups = {}
        for f in glob.glob(os.path.join(monlib_path, "?", "*.cif")):
            doc = gemmi.cif.read(f)
            block = doc.find_block("comp_list")
            for row in block.find("_chem_comp.", ["id", "group"]):
                cgroups[row[0]] = row.str(1)
                try: self.assertTrue(lgroups.get(row[0]) == row.str(1))
                except AssertionError as e: self.errors.append(str(e))

        only_in_cif = set(cgroups) - set(lgroups)
        try: self.assertFalse(only_in_cif, msg="only in cif files")
        except AssertionError as e: self.errors.append(str(e))
        
        only_in_list = set(lgroups) - set(cgroups)
        try: self.assertFalse(only_in_list, msg="only in list")
        except AssertionError as e: self.errors.append(str(e))

    def test_gemmi_monlib(self):
        comp_list, link_list, mod_list = read_mon_lib_list()
        lgroups = {row.str(0):row.str(1) for row in comp_list.find("_chem_comp.", ["id", "group"])}
        try: monlib = gemmi.read_monomer_lib(monlib_path, list(lgroups))
        except Exception as e: self.errors.append(str(e))

    def test_group(self):
        comp_list, link_list, mod_list = read_mon_lib_list()
        lgroups = {row.str(0):row.str(1) for row in comp_list.find("_chem_comp.", ["id", "group"])}
        lg_set = set(lgroups.values())
        all_set = set(["peptide", "P-peptide", "M-peptide", "DNA", "RNA", "pyranose", "ketopyranose", "furanose", "NON-POLYMER"])
        try: self.assertTrue(lg_set.issubset(all_set), msg="groups: {}".format(str(lg_set)))
        except AssertionError as e: self.errors.append(str(e))

    def test_links(self):
        comp_list, link_list, mod_list = read_mon_lib_list()
        lgroups = {row.str(0):row.str(1) for row in comp_list.find("_chem_comp.", ["id", "group"])}
        known_groups = set(lgroups.values())
        known_groups.add("DNA/RNA") # can we really consider these known groups?
        known_groups.add("pept")
        
        ltab = link_list.find("_chem_link.", ["id", "comp_id_1", "mod_id_1", "group_comp_1", "comp_id_2", "mod_id_2", "group_comp_2"])
        mtab = mod_list.find("_chem_mod.", ["id", "comp_id", "group_id"])
        
        link_undef_group = [(row.str(0), row.str(i)) for row in ltab for i in (3,6) if row.str(i) and row.str(i) not in known_groups]
        mod_undef_group =  [(row.str(0), row.str(2)) for row in mtab if row.str(2) and row.str(2) not in known_groups]
        link_undef_comp =  [(row.str(0), row.str(i)) for row in ltab for i in (1,4) if row.str(i) !="" and row.str(i) not in lgroups]
        mod_undef_comp =   [(row.str(0), row.str(1)) for row in mtab if row.str(1) !="" and row.str(1) not in lgroups]

        try: self.assertFalse(link_undef_group, msg="undefined groups referenced in links")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(mod_undef_group, msg="undefined groups referenced in mods")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(link_undef_comp, msg="undefined comp referenced in links")
        except AssertionError as e: self.errors.append(str(e))
        try: self.assertFalse(mod_undef_comp, msg="undefined comp referenced in mods")
        except AssertionError as e: self.errors.append(str(e))

        for row in ltab:
            if row.str(1):
                try: self.assertEqual(row.str(3), lgroups[row.str(1)], msg="{} in link {}".format(row.str(1), row.str(0)))
                except AssertionError as e: self.errors.append(str(e))
            if row.str(4):
                try: self.assertEqual(row.str(6), lgroups[row.str(4)], msg="{} in link {}".format(row.str(4), row.str(0)))
                except AssertionError as e: self.errors.append(str(e))

        for row in mtab:
            if row.str(1):
                try: self.assertEqual(row.str(2), lgroups[row.str(1)], msg="{} in mod {}".format(row.str(1), row.str(0)))
                except AssertionError as e: self.errors.append(str(e))

if __name__ == '__main__':
    unittest.main()
