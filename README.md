# NEDscrape Navigator
A Python tool for retrieving galaxy info from the NASA Extragalactic Database.

Now Python 3 compatible!

Run without any commands for basic usage help::

	python NEDscrape.py

	NEDscrape.py  version: 3.1   04/06/2021

	A program to grab galaxy data from the NED server.

	Normal usage example: python NEDscrape.py -f names.txt -n outFile -o /Users/me/ -m -v

	 -f is a text file of names of objects, one per line, for which data will be retrieved
	 -n is the name of the csv file that will be created and the NED data written to
	 -o is the full pathname describing where the file should be saved
	 -m tells it to use multithreading to speed up data retrieval. Not for large queries!
	 -v tells it to run verbosely
