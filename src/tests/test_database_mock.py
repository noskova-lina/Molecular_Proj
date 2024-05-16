import unittest
from unittest.mock import MagicMock, patch
import sqlite3


class TestDatabaseMock(unittest.TestCase):
    @patch('sqlite3.connect')
    def test_insert_data(self, mock_connect):
        mock_cursor = MagicMock()
        mock_connection = MagicMock()
        mock_connect.return_value = mock_connection
        mock_connection.cursor.return_value = mock_cursor

        conn = sqlite3.connect('dummy')
        cur = conn.cursor()
        cur.execute('INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z) VALUES (?, ?, ?, ?, ?, ?)',
                    ('test_name', 1, 'H', 0.1, 0.2, 0.3))

        mock_cursor.execute.assert_called_with(
            'INSERT INTO molecules (molecule_name, atom_index, atom, x, y, z) VALUES (?, ?, ?, ?, ?, ?)',
            ('test_name', 1, 'H', 0.1, 0.2, 0.3))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestDatabaseMock)
    unittest.TextTestRunner().run(suite)
