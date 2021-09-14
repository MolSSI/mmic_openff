def get_files(ext: str, rel_path: str):
    """Get all files of extension ext in openff data dirs relative to rel_path.
    In the source distribution, these files are in ``openff/toolkit/data/``,
    but on installation, they're moved to somewhere in the user's python
    site-packages directory.

    Parameters
    ----------
    ext : str
        File extension
    rel_path: str
        Path relative to openff data dir: openff/toolkit/data/rel_path

    Returns
    -------
    List of files of extension ext.

    Examples
    --------

    >>> get_files("sdf", "molecules")
    ... [file1.sdf, file2.sdf, ...]
    """

    import os
    import glob
    from pkg_resources import resource_filename

    pathname = resource_filename(
        "openff.toolkit", os.path.join("data", rel_path, "*." + ext)
    )

    return glob.glob(pathname)
