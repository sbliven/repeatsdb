import json
import requests
from bisect import bisect_left
import os
import math
import Bio
import itertools
from Bio import SeqIO
from Bio.PDB import PDBList, MMCIFParser, PDBParser


class RepeatsDBRegion(object): #was PDBRepeats
    """A single repeat from RepeatDB.
    """
    def __init__(self, pdbChain=None, json=None):
        """Initialize a repeat region.
        
        One of pdbChain or json should be given. pdbChain fetches the entry and initializes with the first repeat region.
        The following are equivalent:
        
            RepeatDBRegion(pdbChain)
            RepeatDBRegion(json=RepeatDBRegion.fetchRepeatRegionsJSON(pdbChain)[0])
   
        
        @param pdbChain (str): PDB ID and chain (e.g. '3vszA'). The first repeat region will be used
        @param json (object): JSON specifying a single repeat_region
        """
        if (pdbChain is None) == (json is None): # not exactly one given
            raise Exception("Initialize with either pdbChain or json")
        
        entry = json or self.fetchRepeatRegionsJSON(pdbChain)[0]
        self.json = entry

        # Many elements now one-element lists. Unwrap them
        def unwrap(l):
            assert type(l) == list and len(l) == 1
            return l[0]
        
        self.pfam = entry.get("Pfam",[])
        #self.average_unit = entry["average_unit"]
        self.rdb_class = unwrap(entry["class"])
        self.rdb_classification = unwrap(entry["classification"])
        #self.entry_type = entry["entry_type"]
        #self.has_insertions = entry["has_insertions"]
        #self.has_superegion = entry["has_superegion"]
        self.rdb_id = unwrap(entry["id"])
        self.insertions = entry.get("insertions", [])
        #self.reg_seqres_coverage = entry["reg_seqres_coverage"]
        self.region_id = unwrap(entry["region_id"])
        #self.repDB_source = entry["repDB_source"]
        self.rdb_subclass = unwrap(entry["subclass"])
        self.superregions = entry.get("superreg", [])
        self.units = entry["units"]
        #self.units_number = entry["units_number"]

        self._resseq = None

        # These are essentially tests that I understand RepeatDB data
        assert entry["units_number"] == len(self.units)
        assert len(entry["class"]) == 1
        assert len(entry["subclass"]) == 1
        assert len(entry["classification"]) == 1
        assert entry["has_insertions"] or "insertions" not in entry or len(entry["insertions"]) == 0
        assert entry["has_superegion"] or "superreg" not in entry or len(entry["superreg"]) == 0
        
        # TODO test that this matches. 
        #assert math.floor(float(sum(self.getRepeatLengths())) / len(self.units) + .5) == entry["average_unit"]

    @classmethod
    def fetchRepeatRegionsJSON(cls, pdbChain):
        """Fetch JSON object for all repeat regions in the specified chain
        """
        url = ("http://repeatsdb.bio.unipd.it/ws/search?entry_type=repeat_region"
               "&id={}&collection=repeat_region&show=ALL"
               ).format(pdbChain)

        # Retrive all entries and parse into JSON object
        req = requests.get(url)
        req.raise_for_status()
        entries = json.loads(req.content)

        return entries

    def getRepeatPDBRanges(self):
        """Gets the residues range for each repeat.
        
        Each repeat is defined as an array with alternating [start, end, start, end, ...]
        giving the (inclusive) bounds of each repeat
        """
        ranges = sorted([list(u) for u in self.units])  # deep copy

        # TODO test with insertion codes, missing residues, etc
        i = 0
        for istart, iend in sorted(self.insertions):
            # move i up to range after istart
            while i < len(ranges) and ranges[i][0] <= istart:
                i += 1
            r = ranges[i - 1]
            assert r[0] <= istart and iend <= r[-1]

            # Insert
            j = bisect_left(r, istart)
            r.insert(j, iend+1)
            r.insert(j, istart-1)
        return ranges

    def getRepeatLengths(self):
        lengths = []

        for r in self.getRepeatPDBRanges():
            l = sum([r[i + 1] - r[i] + 1 for i in range(0, len(r), 2)])
            lengths.append(l)
        return lengths

    @property
    def resseq(self):
        """Biopython Seq object representing the RESSEQ for this PDB chain"""
        if not self._resseq:
            pdbId, chain = self.rdb_id[:4], self.rdb_id[4:]
            self._resseq = _fetcher.get_resseq(pdbId, chain)
        return self._resseq


    def _repeat_assignment(self):
        """For each position in the sequence, determines which repeat covers it
        (or None for insertions)

        Stops at the end of the last repeat.

        Returns: generator of int
        """
        ranges = self.getRepeatPDBRanges()

        if not ranges:
            return  # no repeats

        pos = 0  # position in sequence

        for unit in range(len(ranges)):
            for i in range(0, len(ranges[unit]), 2):
                #TODO use SIFT to translate to seqres coordinates
                start = ranges[unit][i]
                end = ranges[unit][i+1]

                # ranges should be sorted
                assert pos <= start
                assert pos <= end

                while pos < start:
                    yield None
                    pos += 1
                while pos <= end:
                    yield unit
                    pos += 1

    def repeatstr(self, width=None):
        seq = str(self.resseq.seq)
        
        # Map repeats to 1-9, A-Z
        chars = [chr(o) for o in itertools.chain(
            range(ord("1"),ord("9")+1),
            range(ord("A"), ord("Z")+1),
            range(ord("a"), ord("z")+1))]

        assignment = "".join([chars[i] if i is not None else " " for i in self._repeat_assignment()])
        assignment += " "*(len(seq)-len(assignment))

        if width is None:
            return "\n".join((seq,assignment))
        
        lines = []
        for i in range(0, len(seq), width):
            lines.append(seq[i:(i + width)])
            lines.append(assignment[i:(i + width)])
        return "\n".join(lines)

    def __len__(self):
        return len(self.units)

    def __str__(self):
        return f"RepeatsDB:{self.rdb_id} with {len(self.units)} repeats"

class Fetcher(object):
    """Fetch structures from the PDB and store them in a RCSB-style directory structure
    
    Compatible with BioJava PDB_DIR environments. 
    """
    def __init__(self, root=None):
        self._pdblists = {}
        self._root = root or os.path.join(os.path.expanduser('~'), 'pdb')

    def get_pdblist(self, file_format="pdb"):
        """Get a PDBList compatible with the RCSB directory structure.
        
        The PDBList can be used to fetch additional structures or access them on disk.
        
        @param file_format (str): type of files to search for (mmCIF|pdb)
        @return Bio.PDBList
        """
        ff = file_format.lower()
        if ff not in self._pdblists:
            if ff == "mmcif":
                prefix = "mmCIF"
            else:
                prefix = ff
            path = os.path.join(self._root, 'pdb', 'data', 'structures', 'divided', prefix)
            self._pdblists[ff] = PDBList(pdb=path)
        return self._pdblists[ff]

    def fetch_structure(self, pdbId, file_format="pdb"):
        """Fetch the specified structure
        
        It will be loaded from disk or downloaded from wwPDB
        """
        filename = self.get_pdblist(file_format).retrieve_pdb_file(pdbId, file_format=file_format)
        if file_format.lower() == "mmcif":
            parser = MMCIFParser()
            struc = parser.get_structure(pdbId, filename)
        elif file_format.lower() == "pdb":
            parser = PDBParser()
            struc = parser.get_structure(pdbId, filename)
        else:
            raise KeyError("Unknown format " + file_format)
        return struc

    def get_resseq(self, pdbId, chain):
        """Parse the seqres for the specified structure and chain
        """
        filename = self.get_pdblist().retrieve_pdb_file(pdbId, file_format="pdb")
        records = [r for r in SeqIO.parse(filename, "pdb-seqres") if r.annotations["chain"] == chain]
        if len(records) < 1:
            raise KeyError("chain %s not found" % chain)
        assert len(records) == 1  # should be one record per chain?
        return records[0]

_fetcher = Fetcher()

class RepeatsDBSearch(object):
    """Interface to RepeatDB Search API
    """

    base_url = "http://repeatsdb.bio.unipd.it"
    def search(self,query,collection="pdb_chain",cluster=None, params=None):
        """
        Low-level access to RepeatsDB search. The databases supports a rich set of query terms,
        which can be combined using boolean logic. 
        
        See http://repeatsdb.bio.unipd.it/help for a full list of supported attributes
        
        Args:
            query (str): query string, e.g. "classification:III.3". Multiple terms may be joined with AND and OR
            collection (str): One of "repeat_region", "pdb_chain", or "uniprot_protein"
            cluster (int): Threshold to filter down to a non-redundant set. Only 40, 60, and 90 are supported.
            params (list of str): list of attributes to include in the result. Default: "ALL"
            
        Return: list of objects representing the JSON response
        """
        # Convert params to comma separated list
        if params is None:
            params = ["ALL"]
        if type(params) == list:
            params = ",".join(params)

        url = ("{base_url}/ws/search?query={query}"
               "&collection={collection}" 
               "&show={params}"
              )
        if cluster:
            if cluster not in (40,60,90):
                raise ValueError("Unsupported clustering threshold. Use 40, 60, 90, or None")
            url += "&filter={cluster}"
        url = url.format( base_url=self.base_url, query=query, collection=collection, params=params, cluster=cluster )
        response = requests.get(url)
        response.raise_for_status()
        # Retrive all entries and parse into JSON object
        entries = json.loads(response.text)
        return entries
    
    def or_search(self,queries,*args,**kwargs):
        """Get union of multiple queries
        
        If needed, queries will be split into multiple requests to fit length limits.
        
        Extra arguments are passed to `search`
        Args:
            queries (list of str): list of query strings
        
        Return: list of objects representing the collective JSON response of each query
        """
        # maximum URL query length per batch
        max_length = 200 #url length
        max_queries = 10 #number of terms
        response = []
        
        curr_len = 0
        curr_queries = []
        for query in queries:
            # Check if there is room for the next query
            if curr_len + len(query) + 2 <= max_length and len(curr_queries) < max_queries or curr_len <= 0:
                curr_queries.append(query)
                curr_len += len(query) + 2
            else:
                # fire request
                response += self.search("OR".join(curr_queries),*args,**kwargs)
                curr_queries = [query]
                curr_len = len(query)
        # fire remaining requests
        if curr_queries:
            response += self.search("OR".join(curr_queries),*args,**kwargs)
            
        #TODO deduplicate?
        return response