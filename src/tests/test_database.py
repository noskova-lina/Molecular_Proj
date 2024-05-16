import unittest
import sqlite3


class TestDatabase(unittest.TestCase):
    def setUp(self):
        self.conn = sqlite3.connect(':memory:')
        self.cur = self.conn.cursor()
        self.cur.execute('''
        CREATE TABLE molecules (
            molecule_name TEXT NOT NULL,
            atom_index INTEGER NOT NULL,
            atom TEXT NOT NULL,
            x REAL NOT NULL,
            y REAL NOT NULL,
            z REAL NOT NULL,
            PRIMARY KEY (molecule_name, atom_index)
        )''')
        self.conn.commit()

    def test_table_created(self):
        self.cur.execute("PRAGMA table_info(molecules)")
        cols = [row[1] for row in self.cur.fetchall()]
        self.assertEqual(cols, ['molecule_name', 'atom_index', 'atom', 'x', 'y', 'z'])

    def test_insert_data(self):
        test_data = ('dsgdb9nsd_000001', 1, 'H', 0.0029, -0.09, 0.0987)
        self.cur.execute('''
        INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z)
        VALUES (?, ?, ?, ?, ?, ?)''', test_data)
        self.conn.commit()

        self.cur.execute('SELECT * FROM molecules WHERE molecule_name = "dsgdb9nsd_000001"')
        self.assertEqual(self.cur.fetchone(), test_data)

    def test_view_data(self):
        test_data = [
            ('dsgdb9nsd_000001', 1, 'H', 0.002150416, -0.0060313176, 0.0019761204),
            ('dsgdb9nsd_000002', 2, 'С', 0.0172574639, -0.00890, 0.0019761204)
        ]
        self.cur.executemany('''
        INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z)
        VALUES (?, ?, ?, ?, ?, ?)''', test_data)
        self.conn.commit()

        self.cur.execute('SELECT * FROM molecules ORDER BY molecule_name')
        fetched_data = self.cur.fetchall()
        self.assertEqual(fetched_data, test_data)

    def test_empty_insert(self):
        with self.assertRaises(sqlite3.IntegrityError):
            self.cur.execute('''
            INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z)
            VALUES (?, ?, ?, ?, ?, ?)''', ('dsgdb9nsd_000003', 1, None, 0.0, 0.0, 0.0))
            self.conn.commit()

    def test_invalid_atom_index(self):
        with self.assertRaises(sqlite3.IntegrityError):
            self.cur.execute('''
            INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z)
            VALUES (?, ?, ?, ?, ?, ?)''', ('dsgdb9nsd_000003', 'invalid', 'H', 0.0, 0.0, 0.0))
            self.conn.commit()

    def test_duplicate_entry(self):
        test_data = ('dsgdb9nsd_000004', 1, 'H', 0.0, 0.0, 0.0)
        self.cur.execute('''
        INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z)
        VALUES (?, ?, ?, ?, ?, ?)''', test_data)
        self.conn.commit()

        with self.assertRaises(sqlite3.IntegrityError):
            self.cur.execute('''
            INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z)
            VALUES (?, ?, ?, ?, ?, ?)''', test_data)
            self.conn.commit()

    def test_fetch_non_existent_data(self):
        self.cur.execute('SELECT * FROM molecules WHERE molecule_name = "non_existent"')
        self.assertIsNone(self.cur.fetchone())

    def tearDown(self):
        self.conn.close()


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDatabase)
    unittest.TextTestRunner().run(suite)
