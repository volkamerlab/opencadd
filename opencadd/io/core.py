"""
opencadd.io.core

Defines the base class for the io module.
"""

from pathlib import Path


class _Base:
    """
    Base class for io classes.
    """

    def __init__(self):
        """
        This class contains no __init__ initialization.
        """

        raise RuntimeError("This class only support initialization from classmethods.")

    @classmethod
    def from_file(cls, filepath, **kwargs):
        """
        Parse a structure from a file into different output data types.

        Parameters
        ----------
        filepath : pathlib.Path or str
            Path to file.
        **kwargs
            Arbitrary keyword arguments.
        """

        raise NotImplementedError("Implement in your subclass!")

    @classmethod
    def from_text(cls, text, **kwargs):
        """
        Parse a structure from a string (text) into different output data types.

        Parameters
        ----------
        text : str
            Structure file text.
        **kwargs
            Arbitrary keyword arguments.
        """

        raise NotImplementedError("Implement in your subclass!")

    @classmethod
    def _file_to_text(cls, filepath):
        """
        Get content (text) from file.

        Parameters
        ----------
        filepath : pathlib.Path or str
            Path to file.

        Returns
        -------
        str
            File content (text).
        """

        filepath = cls._convert_filepath(filepath)

        with open(filepath, "r") as f:
            text = f.read()

        return text

    @classmethod
    def _convert_filepath(cls, filepath):
        """
        Convert a filepath.

        Parameters
        ----------
        filepath : pathlib.Path or str
            Path to file.

        Returns
        -------
        pathlib.Path
            Path to file.

        Raises
        ------
        FileNotFoundError
            Raised if file does not exist.
        """

        filepath = Path(filepath)

        if not filepath.exists():
            raise FileNotFoundError(f"File {filepath} does not exist.")

        return filepath
