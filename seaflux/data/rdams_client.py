#!/usr/bin/env python
"""
Adapted from the rdams_client.py script v2.0.1

written by Doug Schuster (schuster@ucar.edu), Riley Conroy (rpconroy@ucar.edu)
"""
import argparse
import codecs
import getpass
import json
import os
import sys
import pathlib

import requests


BASE_URL = "https://rda.ucar.edu/json_apps/"
USE_NETRC = False
DEFAULT_AUTH_FILE = str(
    pathlib.Path("~/.rdamsrc")
    .expanduser()
    .absolute()
    .resolve()
)


def obfuscate(string):
    """Obfuscate string."""
    return codecs.encode(string, "rot_13")


def unobfuscate(string):
    """Decode obfuscated string."""
    return codecs.decode(string, "rot_13")


def get_userinfo():
    """Get username and password from the command line."""
    user = input("Enter your RDA username or email: ")
    pasw = getpass.getpass("Enter your RDA password: ")
    # try:
    write_pw_file(user, pasw)
    # except Exception as e:
    #    print("Error writing password file: " + str(e))
    return (user, pasw)


def write_pw_file(username, password, pwfile=DEFAULT_AUTH_FILE):
    """Write out file with user information."""
    with open(pwfile, "w") as file_open:
        npwstring = username + "," + password
        ob_str = obfuscate(npwstring)
        file_open.write(ob_str)


def read_pw_file(pwfile):
    """Read user information from pw file.

    Args:
        pwfile (str): location of password file.

    Returns:
        (tuple): (username,password)
    """
    with open(pwfile, "r") as f:
        pwstring = unobfuscate(f.read())
        (username, password) = pwstring.split(",", 2)
    return (username, password)


def check_file_status(filepath, filesize):
    """Prints file download status as percent of file complete.

    Args:
        filepath (str): File being downloaded.
        filesize (int): Expected total size of file in bytes.

    Returns:
        None
    """
    sys.stdout.write("\r")
    sys.stdout.flush()
    size = int(os.stat(filepath).st_size)
    percent_complete = (size / filesize) * 100
    sys.stdout.write("%.3f %s" % (percent_complete, "% Completed"))
    sys.stdout.flush()


def download_files(filelist, out_dir="./", cookie_file=None):
    """Download files in a list.

    Args:
        filelist (list): List of web files to download.
        out_dir (str): directory to put downloaded files

    Returns:
        None
    """
    if cookie_file is None:
        cookies = get_cookies()
    else:
        cookies = cookie_file
    for _file in filelist:
        file_base = os.path.basename(_file)
        out_file = out_dir + file_base
        print("Downloading", file_base)
        req = requests.get(_file, cookies=cookies, allow_redirects=True, stream=True)
        filesize = int(req.headers["Content-length"])
        with open(out_file, "wb") as outfile:
            chunk_size = 1048576
            for chunk in req.iter_content(chunk_size=chunk_size):
                outfile.write(chunk)
                if chunk_size < filesize:
                    check_file_status(out_file, filesize)
        check_file_status(out_file, filesize)
        print()


def get_authentication(pwfile=DEFAULT_AUTH_FILE):
    """Attempts to get authentication.

    Args:
        pwfile (str): location of password file.

    Returns:
        (tuple): username, passord
        (None): If using .netrc file
    """
    if USE_NETRC:
        return None
    if os.path.isfile(pwfile) and os.path.getsize(pwfile) > 0:
        return read_pw_file(pwfile)
    else:
        return get_userinfo()


def get_cookies(username=None, password=None):
    """Authenticates with RDA and returns authentication cookies.

    The user must authenticate with
    authentication cookies per RDA policy.

    Args:
        username (str): RDA username. Typically the user's email.
        password (str): RDA password.

    Returns:
        requests.cookies.RequestsCookieJar: Login request's cookies.
    """
    if username is None and password is None:
        username, password = get_authentication()

    login_url = "https://rda.ucar.edu/cgi-bin/login"
    values = {"email": username, "passwd": password, "action": "login"}
    ret = requests.post(login_url, data=values)
    if ret.status_code != 200:
        print("Bad Authentication")
        print(ret.text)
        exit(1)
    return ret.cookies
