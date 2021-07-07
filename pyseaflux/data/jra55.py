#! /usr/bin/env python
"""
python script to download selected files from rda.ucar.edu
after you save the file, don't forget to make it executable
i.e. - "chmod 755 <name_of_script>"
"""


def download_jra(
    username,
    password,
):
    """Downloads JRA data
    Work in progress"""
    import http.cookiejar as cookiejar
    import os
    import sys

    from urllib import request

    if len(sys.argv) != 2:
        print("usage: " + sys.argv[0] + " [-q] password_on_RDA_webserver")
        print("-q suppresses the progress message for each file that is downloaded")
        sys.exit(1)

    verbose = True
    if len(sys.argv) == 3 and password == "-q":
        verbose = False

    cj = cookiejar.MozillaCookieJar()
    opener = request.build_opener(request.HTTPCookieProcessor(cj))

    # check for existing cookies file and authenticate if necessary
    do_authentication = False
    if os.path.isfile("auth.rda.ucar.edu"):
        cj.load("auth.rda.ucar.edu", False, True)
        for cookie in cj:
            if cookie.name == "sess" and cookie.is_expired():
                do_authentication = True
    else:
        do_authentication = True
    if do_authentication:
        opener.open(
            "https://rda.ucar.edu/cgi-bin/login",
            f"email={username}&password={password}&action=login",
        )
        # save the authentication cookies for future downloads
        # NOTE! - cookies are saved for future sessions because overly-frequent
        # authentication to our server can cause your data access to be blocked
        cj.clear_session_cookies()
        cj.save("auth.rda.ucar.edu", True, True)

    # download the data file(s)
    listoffiles = [
        "anl_surf/1982/anl_surf.001_pres.reg_tl319.1982010100_1982123118",
        "anl_surf/1982/anl_surf.033_ugrd.reg_tl319.1982010100_1982123118",
        "anl_surf/1982/anl_surf.034_vgrd.reg_tl319.1982010100_1982123118",
    ]

    for file in listoffiles:
        idx = file.rfind("/")
        if idx > 0:
            ofile = file[idx + 1 :]
        else:
            ofile = file

        if verbose:
            print("downloading " + ofile + "...")

        infile = opener.open("http://rda.ucar.edu/data/ds628.0/" + file)
        outfile = open(ofile, "wb")

        outfile.write(infile.read())
        outfile.close()

        if verbose:
            print("done")
