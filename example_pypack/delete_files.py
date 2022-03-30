#!/usr/bin/python3

import argparse
import os
import sys
import re

def sanitize_wildcard(wildcard) :
    # Sanitize. That is, any time two unescaped asterisks are next to eachother,
    # fix that.
    wildcard2 = ""
    # Create an array for escaped characters.
    escapes = [False]
    for c in wildcard :
        if c == "\\" :
            escapes.append(not escapes[-1])
        else :
            escapes.append(False)
    # Next, find anytime there are multiple unescaped wildcards.
    asterisks = False
    for i in range(len(wildcard)) :
        if wildcard[i] == "*" and not escapes[i] :
            if asterisks :
                continue
            else :
                asterisks = True
                wildcard2 += "*"
        else :
            asterisks = False
            wildcard2 += wildcard[i]
    return wildcard2

def wildcard_match(string, wildcard) :
    # First, match any non-wildcard things.
    # Base case.
    wildcard2 = sanitize_wildcard(wildcard)
    if len(string) == 0 and len(wildcard2) == 0 :
        return True
    elif len(wildcard2) == 0 :
        return False
    elif len(string) == 0 and any(map(lambda c: c != '*', wildcard2)) :
        return False
    strpos = 0
    wildpos = 0
    escape = False
    while strpos < len(string) and wildpos < len(wildcard2) :
        escape = False
        # If the wildcard has a slash, escape.
        if wildcard2[wildpos] == '\\' :
            escape = True
            wildpos += 1

        # Next, if the wildcard string has a wildcard character,
        # try to match until the final.
        if not escape and wildpos < len(wildcard2) and \
           wildcard2[wildpos] == '*' :
            # Match with the end of the string. We haven't run into errors
            # yet, so the strings must match.
            if wildpos == len(wildcard2) - 1 :
                return True
            elif wildpos < len(wildcard2) - 1 :
                # Go until the string matches the next thing in the
                # wild card.
                while strpos < len(string) and \
                      string[strpos] != wildcard2[wildpos + 1] :
                    strpos += 1
                # Here, the string has not matched the next thing.
                if strpos < len(string) and \
                   string[strpos] != wildcard2[wildpos + 1] :
                    return False
                # From here, split into two.
                res1 = wildcard_match(string[strpos:], wildcard2[wildpos + 1:])
                if strpos < len(string) - 1 :
                    res2 = wildcard_match(string[strpos + 1:],
                                          wildcard2[wildpos:])
                else :
                    res2 = False
                # If one returns true, then the string matches.
                return res1 or res2
        else :
            if string[strpos] != wildcard2[wildpos] :
                return False
            strpos += 1
            wildpos += 1
    if strpos < len(string) or wildpos < len(wildcard2) :
        return False
    return True

def many_wildcards_match(string, wildcards) :
    return any(map(lambda w : wildcard_match(string, w), wildcards))

def find_matches(strings, wildcards) :
    return filter(lambda s : many_wildcards_match(s, wildcards), strings)

def tree(directory = ".") :
    if len(directory) > 0 and directory[-1] == "/" :
        infix = ""
    elif len(directory) > 0 and directory[-1] != "/" :
        infix = "/"
    print(os.listdir(directory))
    entries = [directory + infix + f for f in os.listdir(directory)]
    files = filter(os.path.isfile, entries)
    dirs = filter(os.path.isdir, os.listdir(directory))]
    print(directory)
    print(files)
    print(dirs)
    out = files
    for d in dirs :
        print(d)
        out += [directory + infix + f for f in tree(d)]
    return out

def print_matches(wildcards, recursive = False) :
    files = filter(os.path.isfile, os.listdir())

if __name__ == "__main__" :
    parser = argparse.ArgumentParser(description = """
Deletes files and directories.
""")
    parser.add_argument("--recursive", "-R", action = "store_true",
                        help = """
Recurse through directories to delete items that match within them.
""")
    parser.add_argument("--directories", "-d", action = "store_true",
                        help = """
If directories are specified on the command line, they will have all their
elements deleted, and then will be deleted themselves. Similar to the -r
option for rm.
""")
    parser.add_argument("--print", "-p", action = "store_true",
                        help = """
If this is specified, the names of all the files being deleted will be printed.
""")
    parser.add_argument("files", nargs = '*', help = """
List of files. These may contain wildcard expressions, like *.o, or be a file
name.
""")

    args = vars(parser.parse_args())
    print(os.listdir("./devtools"))
    print(tree())
    

    
