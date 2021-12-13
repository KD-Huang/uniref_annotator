
from utils import say, die, check_path, which, Hit, translate_fasta, try_open
import os

class MMseqs2:
    def __init__(self, c_output_format,
                       c_mmseqs2_filters):

        self.c_output_format = c_output_format
        self.c_mmseqs2_filters = c_mmseqs2_filters


class RunMMseqs2(MMseqs2):

    def __init__(self, c_output_format, c_mmseqs2_filters, mmseqs2 = None, database = None, query = None,
                 seqtype = None, temp = None, mmseqs2_options = None, force_search = None):
        super().__init__(c_output_format, c_mmseqs2_filters)
        self.mmseqs2 = mmseqs2
        self.database = database
        self.query = query
        self.seqtype = seqtype
        self.temp = temp
        self.mmseqs2_options = mmseqs2_options
        self.force_search = force_search
    
    def get_mode(self, path):
        mode = None
        for test in "90 50".split():
            test = "uniref" + test
            if test in path.lower():
                mode = test
        if mode is None:
            die("Could not infer mode from path {}".format(path))

        return mode
    
    def uniref_search(self):

        if which(self.mmseqs2) is None:
            die("<mmseqs2> is not executable as: {}".format(self.mmseqs2))

        for path in [self.database, self.query, self.temp]:
            check_path(path)

        mode = self.get_mode(self.database)
        results = os.path.split(self.query)[1]
        results = os.path.join(self.temp, results)
        results = ".".join([results, mode, "hits"])
        min_id = float(self.get_mode(results).replace("uniref", ""))/100 # converting the minimum id to 1-1.0 range to adapt to mmseqs2
        command = [
            self.mmseqs2,
            "easy-search",
            self.query,
            self.database,
            results,
            self.temp,
            "--format-output", self.c_output_format,
            "--min-seq-id", min_id,
            self.c_mmseqs2_filters
            ]

        command = " ".join(str(k) for k in command)
        command += (" " + self.mmseqs2_options) if self.mmseqs2_options is not None else ""


        if self.force_search or not os.path.exists(results):
            say("Executing:\n ", command)
            os.system(command)
        else:
            say("Using existing results file:\n ", results)

        return results

    def uniref_preidx_search(self, predix):

        if which(self.mmseqs2) is None:
            die("<mmseqs2> is not executable as: {}".format(self.mmseqs2))

        for path in [self.database, self.query, self.temp, predix]:
            check_path(path)

        mode = self.get_mode(self.database)
        results = os.path.split(self.query)[1]
        query_db = os.path.join(self.temp, "qDB."+ mode + "."+ results)
        output_db = os.path.join(self.temp, "resDB." + mode + "." + results)
        results = os.path.join(self.temp, results)
        results = ".".join([results, mode, "hits"])
        min_id = float(self.get_mode(results).replace("uniref", ""))/100
        
        command1 = [self.mmseqs2,
                   "createdb",
                   self.query,
                   query_db
                   ] # prepare a command line for createdb queryDB

        command1 = " ".join(str(k) for k in command1)
        say("Executing:\n ", command1)
        os.system(command1)     
        command2 = [
            self.mmseqs2,
            "search",
            query_db,
            self.database,
            output_db,
            predix,
            "--min-seq-id", min_id,
            self.c_mmseqs2_filters
            ]

        command2 = " ".join(str(k) for k in command2)
        command2 += (" " + self.mmseqs2_options) if self.mmseqs2_options is not None else ""        
        say("Executing:\n ", command2)
        os.system(command2)
        command3 = [self.mmseqs2,
                    "convertalis",
                    query_db,
                    self.database,
                    output_db,
                    results,
                    "--format-output", self.c_output_format,
                    ]
        command3 = " ".join(str(k) for k in command3)
        say("Executing:\n ", command3)
        os.system(command3)

        return results
    def uniref_search_cleanup(self, mmseqs2_hits):
        # this function is to truncate query and subject headers in the mmseqs2 blast-tab output
        # so as to generate consistent sequence identifiers just as those from DIAMOND

        raw_mmseqs2_output = mmseqs2_hits # get default blast tab from mmseqs2
        mmseqs2_blast_tab = open(raw_mmseqs2_output, 'r')
        rowlist = mmseqs2_blast_tab.readlines()
        mmseqs2_blast_tab.close()

        cleaned_mmseqs2_blast_tab = open(raw_mmseqs2_output, 'w')        
        for row in rowlist:
            row_elements = row.rstrip().split('\t')
            cleaned_qseqid = row_elements[0].split()[0] # extract the query identifier excluding explannation
            cleaned_sseqid = row_elements[1].split()[0] # extract the subject identifier excluding the annotation label
            row_elements[0] = cleaned_qseqid # replace old full-length header with the identifier for query
            row_elements[1] = cleaned_sseqid # replace old full-length header with the identifier for query
            cleaned_row = "\t".join(row_elements) + '\n' 
            cleaned_mmseqs2_blast_tab.write(cleaned_row) # overwrite the old result file

        cleaned_mmseqs2_blast_tab.close()

        return raw_mmseqs2_output






        








