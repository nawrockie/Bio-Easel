#include "easel.h"
#include "esl_alphabet.h"
#include "esl_sqio.h"
#include "esl_sq.h"
#include "esl_ssi.h"

/* Macros for converting C structs to perl, and back again)
 * from: http://www.mail-archive.com/inline@perl.org/msg03389.html
 * note the typedef in ~/perl/tw_modules/typedef
 */
#define perl_obj(pointer,class) ({                      \
      SV* ref=newSViv(0); SV* obj=newSVrv(ref, class);  \
      sv_setiv(obj, (IV) pointer); SvREADONLY_on(obj);  \
      ref;                                              \
    })

#define c_obj(sv,type) (                                                \
                        (sv_isobject(sv) && sv_derived_from(sv, #type)) \
                        ? ((type*)SvIV(SvRV(sv)))                       \
                        : NULL                                          \
                                                   )

/* Function:  _c_open_sqfile()
 * Incept:    EPN, Mon Mar  4 13:27:43 2013
 * Synopsis:  Open a sequence file and point a pointer at it.
 * Args:      seqfile:    name of sequence file
 *            do_digital: '1' to read n digital mode, '0' to read in text mode
 *                        digital mode is faster, safer, text preserves case, exact characters
 *            is_rna:     '1' to force RNA alphabet
 *            is_dna:     '1' to force DNA alphabet
 *            is_amino:   '1' to force amino alphabet
 *           
 * Returns:   eslOK on success, dies with 'croak' upon an error.
 */

SV *_c_open_sqfile (char *seqfile, int do_digital, int is_rna, int is_dna, int is_amino)
{
  int           status;      /* Easel status code */
  ESL_SQFILE   *sqfp = NULL; /* open input alignment file */
  /* used only if do_digital is TRUE */
  ESL_ALPHABET *abc       = NULL;
  int           alphatype = eslUNKNOWN;

  /* open input file */
  status = esl_sqfile_Open(seqfile, eslSQFILE_UNKNOWN, NULL, &sqfp);
  if      (status == eslENOTFOUND) croak("Sequence file %s not found.\n",     seqfile);
  else if (status == eslEFORMAT)   croak("Format of file %s unrecognized.\n", seqfile);
  else if (status == eslEINVAL)    croak("Can't autodetect stdin or .gz for sequence file %s\n", seqfile);
  else if (status != eslOK)        croak("Open of sequence file %s failed, code %d.\n", seqfile, status);

  if(do_digital) { 
    if     (is_rna)   alphatype = eslRNA;
    else if(is_dna)   alphatype = eslDNA;
    else if(is_amino) alphatype = eslAMINO;
    else { 
      status = esl_sqfile_GuessAlphabet(sqfp, &alphatype);
      if      (status == eslENOALPHABET) croak("Couldn't guess alphabet from first sequence in %s", seqfile);
      else if (status == eslEFORMAT)     croak("Parse failed (sequence file %s):\n%s\n",
                                               sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));     
      else if (status == eslENODATA)    croak("Sequence file %s contains no data?", seqfile);
      else if (status != eslOK)         croak("Failed to guess alphabet (error code %d)\n", status);
    }
    abc = esl_alphabet_Create(alphatype);
    esl_sqfile_SetDigital(sqfp, abc);
  }

  return perl_obj(sqfp, "ESL_SQFILE");
}    

/* Function:  _c_close_sqfile()
 * Incept:    EPN, Thu Mar 28 10:42:34 2013
 * Synopsis:  Close a sequence file and free the associated ESL_SQFILE.
 * Returns:   eslOK on success, some other status upon failure.
 */

void _c_close_sqfile (ESL_SQFILE *sqfp)
{
  esl_sqfile_Close(sqfp);
  return;
}

/* Function:  _c_open_ssi_index()
 * Incept:    EPN, Mon Mar  4 13:58:29 2013
 * Synopsis:  Open an SSI file for an open sequence file.
 * Returns:   eslOK on success, eslENOTFOUND if SSI file does not exist
 *            dies via croak with informative error message upon an error
 */

int _c_open_ssi_index (ESL_SQFILE *sqfp)
{
  int           status;     /* Easel status code */

  /* Open the SSI index for retrieval */
  if (sqfp->data.ascii.do_gzip)           croak("can't use SSI index for file %s because it is gzipped", sqfp->filename);
  if (esl_sqio_IsAlignment(sqfp->format)) croak("can't use SSI index for file %s because it is an alignment", sqfp->filename);
  
  status = esl_sqfile_OpenSSI(sqfp, NULL);
  
  if      (status == eslEFORMAT)   croak("SSI index for file %s is in incorrect format\n", sqfp->filename);
  else if (status == eslERANGE)    croak("SSI index for file %s is in 64-bit format and we can't read it\n", sqfp->filename);
  else if (status == eslENOTFOUND) return status; /* this is okay, caller may try to deal by creating a new SSI file */
  else if (status != eslOK)        croak("Failed to open SSI index for file %s\n", sqfp->filename);

  return eslOK;
}    

/* Function:  _c_create_ssi_index()
 * Incept:    EPN, Fri Mar  8 09:46:52 2013
 * Synopsis:  Create an SSI index file for an existing sequence file.
 *            Based on and nearly identical to easel's miniapps/esl-sfetch.c::create_ssi_index.
 * Returns:   eslOK on success, eslENOTFOUND if SSI file does not exist
 *            dies via croak with informative error message upon an error
 */

void _c_create_ssi_index (ESL_SQFILE *sqfp)
{
  ESL_NEWSSI *ns      = NULL;
  ESL_SQ     *sq      = NULL; 
  int         nseq    = 0;
  char       *ssifile = NULL;
  uint16_t    fh;
  int         status;

  if(sqfp->do_digital) sq = esl_sq_CreateDigital(sqfp->abc);
  else                 sq = esl_sq_Create();

  esl_strdup(sqfp->filename, -1, &ssifile);
  esl_strcat(&ssifile, -1, ".ssi", 4);
  status = esl_newssi_Open(ssifile, TRUE, &ns); /* TRUE is for allowing overwrite. */
  if      (status == eslENOTFOUND)   croak("failed to open SSI index %s", ssifile);
  else if (status == eslEOVERWRITE)  croak("SSI index %s already exists; delete or rename it", ssifile); /* won't happen, see TRUE above... */
  else if (status != eslOK)          croak("failed to create a new SSI index");

  if (esl_newssi_AddFile(ns, sqfp->filename, sqfp->format, &fh) != eslOK)
    croak("Failed to add sequence file %s to new SSI index\n", sqfp->filename);

  while ((status = esl_sqio_ReadInfo(sqfp, sq)) == eslOK)
    {
      nseq++;
      if (sq->name == NULL) croak("Every sequence must have a name to be indexed. Failed to find name of seq #%d\n", nseq);

      if (esl_newssi_AddKey(ns, sq->name, fh, sq->roff, sq->doff, sq->L) != eslOK)
	croak("Failed to add key %s to SSI index", sq->name);

      if (sq->acc[0] != '\0') {
	if (esl_newssi_AddAlias(ns, sq->acc, sq->name) != eslOK)
	  croak("Failed to add secondary key %s to SSI index", sq->acc);
      }
      esl_sq_Reuse(sq);
    }
  if      (status == eslEFORMAT) croak("Parse failed (sequence file %s):\n%s\n",
					   sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (status != eslEOF)     croak("Unexpected error %d reading sequence file %s",
					    status, sqfp->filename);

  /* Determine if the file was suitable for fast subseq lookup. */
  if (sqfp->data.ascii.bpl > 0 && sqfp->data.ascii.rpl > 0) {
    if ((status = esl_newssi_SetSubseq(ns, fh, sqfp->data.ascii.bpl, sqfp->data.ascii.rpl)) != eslOK) 
      croak("Failed to set %s for fast subseq lookup.", sqfp->filename);
  }

  /* Save the SSI file to disk */
  if (esl_newssi_Write(ns) != eslOK)  croak("Failed to write keys to ssi file %s\n", ssifile);

  /* done */
  esl_sqfile_Position(sqfp, 0); /* rewind b/c we're at the end of the file, and if we try to read it we'll get EOF */
  free(ssifile);
  esl_sq_Destroy(sq);
  esl_newssi_Close(ns);
}    

/* Function:  _c_fetch_one_sequence()
 * Incept:    EPN, Tue Sep 25 17:29:17 2018
 * Purpose:   Fetch a single sequence named <sqname> from 
 *            <sqfp> into <sq>.
 *            As a special case, if <sqname> is NULL, 
 *            we fetch the next sequence in the file.
 *
 * Args:      sqfp   - open ESL_SQFILE to fetch seq from
 *            sqname - name of sequence to fetch
 *            ret_sq - ESL_SQ object to fetch sequence into, created here
 *
 * Returns:   void
 *
 * Dies:      with croak if 
 *            - sequence <sqname> does not exist in <sqfp>
 *            - problem parsing SSI index for <sqfp>
 *            - unable to fetch sequence for some reason
 *            - we fetch the wrong sequence for some reason
 */
void _c_fetch_one_sequence(ESL_SQFILE *sqfp, char *sqname, ESL_SQ **ret_sq) { 
  int      status;    /* Easel status code */
  ESL_SQ  *sq = NULL; /* the sequence */

  /* make sure SSI is valid if we're going to use it */
  if ((sqname != NULL) && (sqfp->data.ascii.ssi == NULL)) croak("sequence file has no SSI information\n", sqfp->filename); 

  if(sqfp->do_digital) sq = esl_sq_CreateDigital(sqfp->abc);
  else                 sq = esl_sq_Create();

  /* from esl-sfetch.c's onefetch() */
  if(sqname != NULL) { 
    status = esl_sqfile_PositionByKey(sqfp, sqname);
    if      (status == eslENOTFOUND) croak("seq %s not found in SSI index for file %s\n", sqname, sqfp->filename); 
    else if (status == eslEFORMAT)   croak("Failed to parse SSI index for %s\n", sqfp->filename);
    else if (status != eslOK)        croak("Failed to look up location of seq %s in SSI index of file %s\n", sqname, sqfp->filename);
  }
  /* if sqname is NULL, we just read the next sequence */

  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) croak("Parse failed (sequence file %s):\n%s\n",  sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (status == eslEOF)     croak("Unexpected EOF reading sequence file %s\n", sqfp->filename);
  else if (status != eslOK)      croak("Unexpected error %d reading sequence file %s\n", status, sqfp->filename);

  if (sqname != NULL && strcmp(sqname, sq->name) != 0 && strcmp(sqname, sq->acc) != 0) 
    croak("whoa, internal error; found the wrong sequence %s, not %s\n", sq->name, sqname);

  *ret_sq = sq;

  return;
}

/* Function:  _c_fetch_one_subsequence()
 * Incept:    EPN, Tue Sep 25 18:30:28 2018
 * Purpose:   Fetch a single subsequence from <start>..<end>
 *            for a sequence named <sqname> from 
 *            <sqfp> into <sq>.
 *            As a special case, if <sqname> is NULL, 
 *            we fetch the next sequence in the file.
 *
 * Args:      sqfp   - open ESL_SQFILE to fetch seq from
 *            sqname - name of sequence to fetch
 *            newname      - name to assign to the sequence, if null <sqname><given_start>-<given_end> will be used
 *            given_start - start position passed in (may be > given_end, if we want the revcomp)
 *            given_end   - end position passed in (may be < given_start, if we want the revcomp)
 *            do_res_revcomp - TRUE to force revcomp of a length 1 sequence, since
 *                             it's impossible to tell from given_start/given_end
 *                             if a 1 residue sequence should be revcomp'ed.
 *            ret_sq - ESL_SQ object to fetch sequence into, created here
 *
 * Returns:   void
 *
 * Dies:      with croak if 
 *            - sequence <sqname> does not exist in <sqfp>
 *            - problem parsing SSI index for <sqfp>
 *            - unable to fetch subsequence for some reason
 */
void _c_fetch_one_subsequence(ESL_SQFILE *sqfp, char *sqname, char *newname, long given_start, long given_end, int do_res_revcomp, ESL_SQ **ret_sq) { 
  ESL_SQ  *sq = NULL; /* the sequence */
  int     start, end;            /* start/end for esl_sqio_FetchSubseq() */
  int     do_revcomp;            /* are we revcomp'ing? */

  /* make sure SSI is valid if we're going to use it */
  if (sqfp->data.ascii.ssi == NULL) croak("sequence file has no SSI information\n"); 

  if(sqfp->do_digital) sq = esl_sq_CreateDigital(sqfp->abc);
  else                 sq = esl_sq_Create();

  /* reverse complement indicated by coords. */
  if (given_end != 0 && given_start > given_end) { 
    start = given_end; end   = given_start; do_revcomp = TRUE; 
  }
  else if (given_end == given_start && do_res_revcomp) { 
    /* odd case: single residue, can't tell from given_start/given_end if we should revcomp */
    start = given_end; end = given_start;   do_revcomp = TRUE;  
  }
  else { 
    start = given_start; end = given_end;   do_revcomp = FALSE; 
  }

  /* fetch the subsequence, croak upon an error */
  if (esl_sqio_FetchSubseq(sqfp, sqname, start, end, sq) != eslOK) croak(esl_sqfile_GetErrorBuf(sqfp));

  if      (newname != NULL) esl_sq_SetName(sq, newname);
  else                      esl_sq_FormatName(sq, "%s/%d-%d", sqname, given_start, (given_end == 0) ? sq->L : given_end);

  /* possibly reverse complement the subseq we just fetched */
  if (do_revcomp) { 
    if (esl_sq_ReverseComplement(sq) != eslOK) croak("Failed to reverse complement %s; is it a protein?\n", sq->name);
  }

  *ret_sq = sq;

  return;
}

/* Function:  _c_sq_to_seqstring()
 * Incept:    EPN, Mon Mar 25 15:35:13 2013
 * Synopsis:  Construct a sequence string from an ESL_SQ.
 * Args:      sq    - the ESL_SQ object
 *            textw - width for each sequence of FASTA record, -1 for unlimited.
 *            key   - key used to fetch sequence by caller, useful only for informative error output
 * Returns:   A pointer to a string that is the sequence in FASTA format.
 * Dies:      With croak, if we run out of memory or if the sequence is digitized
 *            and we are unable to textize it.
 */
char *_c_sq_to_seqstring (ESL_SQ *sq, int textw, char *key, int64_t *ret_n)
{    
  int     status = eslOK;        /* easel return status */
  char   *seqstring = NULL;      /* the sequence string */
  int64_t n   = 0;               /* position in string */
  int64_t n2  = 0;               /* position in string */
  int64_t pos = 0;               /* position in sq->seq */

  if(sq->dsq) { /* sequence is digitized, textize it */
    if((status = esl_sq_Textize(sq)) != eslOK) croak("problem converting digitized sequence to text sequence (%s)\n", key);
  }

  /* '>' character */
  n = 0;
  if (esl_strdup(">", 1, &seqstring)              != eslOK) croak("out of memory while fetching sequence %s\n", key);
  n++;
  /* name */
  n2 = strlen(sq->name);
  if (esl_strcat(&seqstring, n, sq->name, n2)     != eslOK) croak("out of memory while fetching sequence %s\n", key);
  n += n2;
  /* accession */
  if (sq->acc[0] != 0) { 
    if (esl_strcat(&seqstring, n, " ", 1)         != eslOK) croak("out of memory while fetching sequence %s\n", key);
    n++;
    n2 = strlen(sq->acc);
    if (esl_strcat(&seqstring, n, sq->acc, n2)    != eslOK) croak("out of memory while fetching sequence %s\n", key);
    n += n2;
  }
  /* desc */
  if (sq->desc[0] != 0) { 
    if (esl_strcat(&seqstring, n, " ", 1)         != eslOK) croak("out of memory while fetching sequence %s\n", key);
    n++;
    n2 = strlen(sq->desc);
    if (esl_strcat(&seqstring, n, sq->desc, n2)   != eslOK) croak("out of memory while fetching sequence %s\n", key);
    n += n2;
  }
  /* newline */
  if (esl_strcat(&seqstring, n, "\n", 1)          != eslOK) croak("out of memory while fetching sequence %s\n", key);
  n++;
  /* sequence */
  if (textw == -1) { /* unlimited line length */
    if (esl_strcat(&seqstring, n, sq->seq, sq->n) != eslOK) croak("out of memory while fetching sequence %s\n", key);
    n += sq->n;
    /* newline */
    if (esl_strcat(&seqstring, n, "\n", 1)        != eslOK) croak("out of memory while fetching sequence %s\n", key);
    n++;
  }
  else { /* textw != -1, limit line length to textw */
    for (pos = 0; pos < sq->n; pos += textw) { 
      n2 = ESL_MIN(textw, sq->n - pos);
      if (esl_strcat(&seqstring, n, sq->seq+pos, n2) != eslOK) croak("out of memory while fetching sequence %s\n", key);
      n += n2;
      /* newline */
      if (esl_strcat(&seqstring, n, "\n", 1) != eslOK) croak("out of memory while fetching sequence %s\n", key);
      n++;
    }
  }

  if(ret_n != NULL) *ret_n = n;
  
  return seqstring;
}

/* Function:  _c_fetch_seq_to_fasta_string()
 * Incept:    EPN, Mon Mar  4 15:03:11 2013
 * Synopsis:  Fetch a sequence from an open sequence file and return it as a FASTA
 *            formatted string.
 * Args:      sqfp  - open ESL_SQFILE to fetch seq from
 *            key   - name or accession of sequence to fetch
 *            textw - width for each sequence of FASTA record, -1 for unlimited.
 * Returns:   A pointer to a string that is the sequence in FASTA format.
 * Dies:      if problem reading sequence
 */
SV *_c_fetch_seq_to_fasta_string (ESL_SQFILE *sqfp, char *key, int textw)
{
  ESL_SQ *sq = NULL;             /* the sequence */
  char   *seqstring = NULL;      /* the sequence string */
  SV     *seqstringSV;           /* SV version of seqstring */
  int64_t n;                     /* length of seqstring */

  /* make sure textw makes sense and we're not in digital mode */
  if(textw < 0 && textw != -1) croak("invalid value for textw\n"); 

  /* fetch the sequence, this will die with croak if there's a problem */
  _c_fetch_one_sequence(sqfp, key, &sq);

  seqstring = _c_sq_to_seqstring(sq, textw, key, &n);
  esl_sq_Destroy(sq);

  seqstringSV = newSVpv(seqstring, n);
  free(seqstring);

  return seqstringSV;
}

/* Function:  _c_fetch_seq_to_fasta_string_given_ssi_number()
 * Incept:    EPN, Mon Apr 23 11:18:17 2018
 * Synopsis:  Fetch a sequence from an open sequence file given its
 *            SSI index number (see _c_fetch_seq_to_fasta_string() for
 *            fetching by key (seqname or accn)) and return it as a FASTA
 *            formatted string.
 * Args:      sqfp  - open ESL_SQFILE to fetch seq from
 *            nkey  - index of the key to return
 *            textw - width for each sequence of FASTA record, -1 for unlimited.
 * Returns:   A pointer to a string that is the sequence in FASTA format.
 * Dies:      if problem reading sequence
 */
SV *_c_fetch_seq_to_fasta_string_given_ssi_number (ESL_SQFILE *sqfp, int nkey, int textw)
{

  int     status;                /* Easel status code */
  ESL_SQ *sq = NULL;             /* the sequence */
  char   *seqstring = NULL;      /* the sequence string */
  SV     *seqstringSV;           /* SV version of seqstring */
  int64_t n;                     /* length of seqstring */

  if(sqfp->do_digital) sq = esl_sq_CreateDigital(sqfp->abc);
  else                 sq = esl_sq_Create();

  /* make sure textw makes sense and we're not in digital mode */
  if(textw < 0 && textw != -1) croak("invalid value for textw\n"); 

  /* adapted from esl-sfetch.c's onefetch() for PositionByNumber instead of PositionByKey */
  if (sqfp->data.ascii.ssi == NULL) croak("sequence file has no SSI information\n"); 
  status = esl_sqfile_PositionByNumber(sqfp, nkey);
  if      (status == eslENOTFOUND) croak("seq index %d not found in SSI index for file %s\n", nkey, sqfp->filename); 
  else if (status == eslEFORMAT)   croak("Failed to parse SSI index for %s\n", sqfp->filename);
  else if (status != eslOK)        croak("Failed to look up location of seq index %d in SSI index of file %s\n", nkey, sqfp->filename);

  status = esl_sqio_Read(sqfp, sq);
  if      (status == eslEFORMAT) croak("Parse failed (sequence file %s):\n%s\n",  sqfp->filename, esl_sqfile_GetErrorBuf(sqfp));
  else if (status == eslEOF)     croak("Unexpected EOF reading sequence file %s\n", sqfp->filename);
  else if (status != eslOK)      croak("Unexpected error %d reading sequence file %s\n", status, sqfp->filename);

  seqstring = _c_sq_to_seqstring(sq, textw, sq->name, &n);
  esl_sq_Destroy(sq);

  seqstringSV = newSVpv(seqstring, n);
  free(seqstring);

  return seqstringSV;
}

/* Function:  _c_fetch_next_seq_to_fasta_string()
 * Incept:    EPN, Fri Mar  8 05:51:49 2013
 * Synopsis:  Fetch the next sequence from an open sequence file and return it as a 
 *            FASTA formatted string. This is a wrapper for _c_fetch_seq_to_fasta_string
 *            that passes NULL for <key>. This function is only necessary because
 *            I don't know how to pass a NULL value in for a char * to an inline C 
 *            function from Perl (as far as I can tell, it can't be done).
 * Args:      sqfp  - open ESL_SQFILE to fetch seq from
 *            textw - width for each sequence of FASTA record, -1 for unlimited.
 * Returns:   A pointer to a string that is the sequence in FASTA format.
 */

SV *_c_fetch_next_seq_to_fasta_string (ESL_SQFILE *sqfp, int textw)
{
  return _c_fetch_seq_to_fasta_string(sqfp, NULL, textw);
}

/* Function:  _c_fetch_subseq_to_fasta_string()
 * Incept:    EPN, Sat Mar 23 05:34:15 2013
 * Synopsis:  Fetch a subsequence.
 *
 * Purpose:   Fetch a subsequence from an open sequence file and return
 *            it as a FASTA formatted string. Based on esl-sfetch's
 *            onefetch_subseq(). The subsequence fetched is from
 *            position <given_start> to <given_end>, as a special case
 *            if <given_end> is 0, then the subsequence is fetched all
 *            the way to the end. If <given_start> > <given_end> (and
 *            <given_end> != 0) the caller is indicating they want the
 *            reverse complement of the subsequence from <given_end>
 *            to <given_start>, we'll fetch the top strand
 *            subsequence, then revcomp it, then return it.
 *
 * Args:      sqfp        - open ESL_SQFILE to fetch seq from
 *            key         - name or accession of sequence to fetch
 *            newname     - name to assign to fetched subsequence
 *            given_start - first position of subseq
 *            given_end   - final position of subseq 
 *            textw - width for each sequence of FASTA record, -1 for unlimited.
 *            do_res_revcomp - TRUE to force revcomp of a length 1 sequence, since
 *                             it's impossible to tell from given_start/given_end
 *                             if a 1 residue sequence should be revcomp'ed.
 *
 * Returns:   A pointer to a string that is the subsequence in FASTA format.
 */

SV *_c_fetch_subseq_to_fasta_string (ESL_SQFILE *sqfp, char *key, char *newname, long given_start, long given_end, int textw, int do_res_revcomp)
{
  ESL_SQ *sq = NULL;             /* the sequence */
  char   *seqstring = NULL;      /* the sequence string */
  SV     *seqstringSV;           /* SV version of seqstring */
  int64_t n;                     /* length of seqstring */

  /* make sure textw makes sense */
  if(textw < 0 && textw != -1) croak("invalid value for textw\n"); 
  /* make sure SSI is valid */
  if (sqfp->data.ascii.ssi == NULL) croak("sequence file has no SSI information\n"); 

  _c_fetch_one_subsequence(sqfp, key, newname, given_start, given_end, do_res_revcomp, &sq);

  seqstring = _c_sq_to_seqstring(sq, textw, key, &n);
  esl_sq_Destroy(sq);

  seqstringSV = newSVpv(seqstring, n);
  free(seqstring);

  return seqstringSV;

}
/* Function:  _c_fetch_seq_name_and_length_given_ssi_number()
 * Incept:    EPN, Mon Apr  8 09:34:01 2013
 * Purpose:   Fetch the primary key of a sequence and the sequence length, 
 *            given it's rank in the SSIindex for the file. This will not 
 *            return the name and length of the <nkey>'th sequence in the 
 *            file, but rather the <nkey>'th sequence as it's ordered in 
 *            the SSI file. 
 *
 *            The return string is a concatenation of the primary key and
 *            the sequence length, separated by a single space (' '). This
 *            is done so we don't have to return two separate values.
 *            
 * Args:      sqfp        - open ESL_SQFILE to fetch seq from
 *            nkey        - index of the key to return
 *
 * Returns:   A pointer to a string that is primary key of the <nkey>'th key
 *            followed by a single ' ', followed by the length of the sequence.
 *            Note that the key cannot contain any spaces, so the only ' '
 *            in the returned string will be the divider between the key
 *            and the sequence length.
 */

SV *_c_fetch_seq_name_and_length_given_ssi_number(ESL_SQFILE *sqfp, int nkey) { 
  int     status;             /* Easel status code */
  char   *key_and_L   = NULL; /* key and L string to return */
  SV     *key_and_LSV;        /* SV version of key and L string to return */
  int64_t L;                  /* length of sequence */
  int64_t Ldup;               /* copy of L we can modify */
  int     Lwidth = 0;         /* number of characters we need for a string conversion of L */
  char   *Lstr = NULL;        /* string of just L */

  /* make sure SSI is valid */
  if (sqfp->data.ascii.ssi == NULL) croak("sequence file has no SSI information\n"); 

  /* fetch the seq name */
  status = esl_ssi_FindNumber(sqfp->data.ascii.ssi, nkey, NULL, NULL, NULL, &L, &key_and_L);
  if     (status == eslEMEM)      croak("out of memory");
  else if(status == eslENOTFOUND) croak("there is no sequence %d\n", nkey);
  else if(status == eslEFORMAT)   croak("error fetching sequence num %d, something wrong with SSI index?\n", nkey);
  else if(status != eslOK)        croak("error fetching sequence num %d\n", nkey);

  Ldup = L;
  Lwidth = 1;
  while(Ldup >= 10) { Ldup/=10; Lwidth++; }
  ESL_ALLOC(Lstr, sizeof(char) * (Lwidth + 1));
  snprintf(Lstr, Lwidth+1, "%" PRId64 "", L);
  
  /* add " " . L to end of key_and_L */
  if((status = esl_strcat(&key_and_L, -1, " ", 1))        != eslOK) croak("out of memory");
  if((status = esl_strcat(&key_and_L, -1, Lstr,  Lwidth)) != eslOK) croak("out of memory");
  free(Lstr);

  key_and_LSV = newSVpv(key_and_L, strlen(key_and_L));
  free(key_and_L);

  return key_and_LSV;

 ERROR: 
  croak("out of memory");
  return NULL; /* NEVER REACHED */
}

/* Function:  _c_fetch_seq_length_given_name()
 * Incept:    EPN, Mon Nov 25 05:09:35 2013
 * Purpose:   Fetch the length of a sequence given its name (primary key).
 *
 *            If the fetched length is 0, then the lengths are unset in
 *            the SSI file. Caller must deal with this.
 *            
 * Args:      sqfp   - open ESL_SQFILE to fetch seq from
 *            sqname - name of sequence we want the length of
 *
 * Returns:   the length of the sequence named <sqname> in <sqfp>,
 *            '0' if lengths are unset in <sqfp>
 *            '-1' if no sequence named <sqname> exists in <sqfp>
 * Dies:      if out of memory or somethings wrong with the SSI file
 */

long _c_fetch_seq_length_given_name(ESL_SQFILE *sqfp, char *sqname) { 
  int      status;   /* Easel status code */
  int64_t  L;        /* length of sequence */
  uint16_t fh;       /* file handle sequence is in, irrelevant since we only have 1 file */
  off_t    roff;     /* offset of start of sqname's record, irrelevant here */

  /* make sure SSI is valid */
  if (sqfp->data.ascii.ssi == NULL) croak("sequence file has no SSI information\n"); 

  /* fetch the length */
  status = esl_ssi_FindName(sqfp->data.ascii.ssi, sqname, &fh, &roff, NULL, &L);
  if     (status == eslEMEM)      croak("out of memory");
  else if(status == eslENOTFOUND) return -1; /* if seq does not exist, we don't die but return -1, caller must deal */
  else if(status == eslEFORMAT)   croak("error fetching sequence name %s, something wrong with SSI index?\n", sqname);
  else if(status != eslOK)        croak("error fetching sequence name %s\n", sqname);

  return L;
}

/* Function:  _c_check_seq_exists()
 * Incept:    EPN, Wed Sep 19 11:34:12 2018
 * Purpose:   Check if a sequence exists given its name (primary key).
 * Args:      sqfp   - open ESL_SQFILE to fetch seq from
 *            sqname - name of sequence we want the length of
 *
 * Returns:   '1' if sequence <sqname> exists in <sqfp>,
 *            '0' if sequence <sqname> does not exist in <sqfp>
 * Dies:      if out of memory or somethings wrong with the SSI file
 */

int _c_check_seq_exists(ESL_SQFILE *sqfp, char *sqname) { 
  int      status;   /* Easel status code */
  uint16_t fh;       /* file handle sequence is in, irrelevant since we only have 1 file */
  off_t    roff;     /* offset of start of sqname's record, irrelevant here */

  /* make sure SSI is valid */
  if (sqfp->data.ascii.ssi == NULL) croak("sequence file has no SSI information\n"); 

  /* look-up sequence */
  status = esl_ssi_FindName(sqfp->data.ascii.ssi, sqname, &fh, &roff, NULL, NULL);
  if     (status == eslEMEM)      croak("out of memory");
  else if(status == eslEFORMAT)   croak("error fetching sequence name %s, something wrong with SSI index?\n", sqname);
  else if(status == eslENOTFOUND) return 0; /* seq does not exist, return '0' */
  else if(status != eslOK)        croak("error fetching sequence name %s\n", sqname);

  return 1; /* if we get here, seq exists */
}

/* Function:  _c_compare_seq_to_seq()
 * Incept:    EPN, Tue Sep 25 10:51:38 2018
 * Purpose:   Check if a sequence exists and is identical
 *            (in residues only) in two sequence files.
 * Args:      sqfp1   - first  open ESL_SQFILE to fetch seq from
 *            sqfp2   - second open ESL_SQFILE to fetch seq from
 *            sqname1 - name of sequence in sqfp1 we want to compare
 *            sqname2 - name of sequence in sqfp2 we want to compare
 *
 * Returns:   '1' if sequence <sqname1> exists in <sqfp1> and 
 *            <sqname2> exists in <sqfp2> and the seqs are 
 *            identical.
 *            '0' if sequence <sqname> exists in <sqfp1> and 
 *            <sqname2> exists in <sqfp2> and the seqs are 
 *            not identical.
 * Dies:      - if sequence <sqname> does not exist in either <sqfp1>
 *              of <sqfp2>
 *            - if we run out of memory
 *            - if both files are not digitized or both not digitized
 *            - something is wrong with the SSI files
 */

int _c_compare_seq_to_seq(ESL_SQFILE *sqfp1, ESL_SQFILE *sqfp2, char *sqname1, char *sqname2) { 
  ESL_SQ     *sq1 = NULL; /* sequence read from first sequence file */
  ESL_SQ     *sq2 = NULL; /* sequence read from second sequence file */

  /* make sure SSI is valid */
  if (sqfp1->data.ascii.ssi == NULL) croak("sequence file 1 %s has no SSI information\n", sqfp1->filename); 
  if (sqfp2->data.ascii.ssi == NULL) croak("sequence file 2 %s has no SSI information\n", sqfp2->filename); 

  /* make sure both sequence files are either both in digital mode or both in text mode */
  if (sqfp1->do_digital && (! sqfp2->do_digital)) croak("sequence file 1 %s is digitized, but sequence file 2 is not digitized %s\n", sqfp1->filename, sqfp2->filename); 
  if (sqfp2->do_digital && (! sqfp1->do_digital)) croak("sequence file 2 %s is digitized, but sequence file 1 is not digitized %s\n", sqfp2->filename, sqfp1->filename); 

  /* if we get here, either both sequence files are digitized or both are not */
  _c_fetch_one_sequence(sqfp1, sqname1, &sq1);
  _c_fetch_one_sequence(sqfp2, sqname2, &sq2);

  /* if the sequences are not the same length, they can't be equal, return 0 now */
  if (sq1->n != sq2->n) return 0;

  /* compare sequences, either digitized or text */
  if(sq1->dsq && sq2->dsq) {
    if (memcmp(sq1->dsq, sq2->dsq, sizeof(ESL_DSQ) * (sq1->n+2)) != 0) return 0; /* return 0 if they're not identical */
  }
  else if (sq1->seq && sq2->seq) { /* sequences are in text mode */
    if (strcmp(sq1->seq, sq2->seq) != 0) return 0; /* return 0 if they're not identical */
  }
  else { 
    croak("whoa, internal error, sequence file types matched but both seqs are not dsq and both seqs are not text\n");
  }

  return 1; /* if we get here, sequences have been compared and are identical */
}

/* Function:  _c_compare_seq_to_subseq()
 * Incept:    EPN, Tue Sep 25 18:22:41 2018
 * Purpose:   Check if a sequence and subsequence exist
 *            and are identical (in residues only) in two sequence files.
 * 
 * Args:      sqfp1   - first  open ESL_SQFILE to fetch seq <sqname1> from
 *            sqfp2   - second open ESL_SQFILE to fetch seq <sqname2> from
 *            sqname1 - name of sequence we want from sqfp1
 *            sqname2 - name of sequence we want from sqfp2
 *            start2  - start position of subsequence 2
 *            end2    - end position of subsequence 2
 * 
 * Returns:   '1' if sequence <sqname> exists in <sqfp1> and <sqfp2>
 *            and is identical.
 *            '0' if sequence <sqname> exists in <sqfp1> and <sqfp2>
 *            and is not identical.
 * Dies:      - if sequence <sqname> does not exist in either <sqfp1>
 *              of <sqfp2>
 *            - if we run out of memory
 *            - if both files are not digitized or both not digitized
 *            - something is wrong with the SSI files
 */

int _c_compare_seq_to_subseq(ESL_SQFILE *sqfp1, ESL_SQFILE *sqfp2, char *sqname1, char *sqname2, long start2, long end2) { 
  ESL_SQ     *sq1 = NULL; /* sequence read from first sequence file */
  ESL_SQ     *sq2 = NULL; /* sequence read from second sequence file */

  /* make sure SSI is valid */
  if (sqfp1->data.ascii.ssi == NULL) croak("sequence file 1 %s has no SSI information\n", sqfp1->filename); 
  if (sqfp2->data.ascii.ssi == NULL) croak("sequence file 2 %s has no SSI information\n", sqfp2->filename); 

  /* make sure both sequence files are either both in digital mode or both in text mode */
  if (sqfp1->do_digital && (! sqfp2->do_digital)) croak("sequence file 1 %s is digitized, but sequence file 2 is not digitized %s\n", sqfp1->filename, sqfp2->filename); 
  if (sqfp2->do_digital && (! sqfp1->do_digital)) croak("sequence file 2 %s is digitized, but sequence file 1 is not digitized %s\n", sqfp2->filename, sqfp1->filename); 

  /* if we get here, either both sequence files are digitized or both are not */
  _c_fetch_one_sequence(sqfp1, sqname1, &sq1);
  _c_fetch_one_subsequence(sqfp2, sqname2, NULL, start2, end2, /*do_res_revcomp=*/0, &sq2);

  /* if the sequences are not the same length, they can't be equal, return 0 now */
  if (sq1->n != sq2->n) return 0;

  /* compare sequences, either digitized or text */
  if(sq1->dsq && sq2->dsq) {
    if (memcmp(sq1->dsq, sq2->dsq, sizeof(ESL_DSQ) * (sq1->n+2)) != 0) return 0; /* return 0 if they're not identical */
  }
  else if (sq1->seq && sq2->seq) { /* sequences are in text mode */
    if (strcmp(sq1->seq, sq2->seq) != 0) return 0; /* return 0 if they're not identical */
  }
  else { 
    croak("whoa, internal error, sequence file types matched but both seqs are not dsq and both seqs are not text\n");
  }

  return 1; /* if we get here, sequences have been compared and are identical */
}


/* Function:  _c_nseq_ssi
 * Incept:    EPN, Mon Apr  8 13:05:39 2013
 * Purpose:   Return the number of sequences in a sequence file.
 *            
 * Args:      sqfp - open ESL_SQFILE
 *
 * Returns:   Number of sequences in the file.
 */

long _c_nseq_ssi(ESL_SQFILE *sqfp) { 
  /* make sure SSI is valid */
  if (sqfp->data.ascii.ssi == NULL) croak("sequence file has no SSI information\n"); 

  return sqfp->data.ascii.ssi->nprimary;
}

/* Function:  _c_nres_ssi
 * Incept:    EPN, Thu Dec  5 06:14:02 2013
 * Purpose:   Return the number of residues in a sequence file.
 *            
 * Args:      sqfp - open ESL_SQFILE
 *
 * Returns:   Number of residues in the file.
 */

SV * _c_nres_ssi(ESL_SQFILE *sqfp) { 
  int status;
  int i;                    /* counter over seqs */
  int64_t nseq;             /* number of seqs in sqfp */
  int64_t nres = 0;         /* summed length of all seqs */
  int64_t nres_dup;         /* duplicate of nres, used to create nres_str */
  char   *nres_str = NULL;  /* string version of nres, to return */
  int     nres_wid = 0;     /* num digits in nres */
  SV     *nres_str_SV;      /* SV version of nres_str, to return */
  int64_t L;                /* length of current seq */

  /* make sure SSI is valid */
  if (sqfp->data.ascii.ssi == NULL) croak("sequence file has no SSI information\n"); 
  
  nseq = sqfp->data.ascii.ssi->nprimary;
  
  for(i = 0; i < nseq; i++) { 
    status = esl_ssi_FindNumber(sqfp->data.ascii.ssi, i, NULL, NULL, NULL, &L, NULL);
    if     (status == eslEMEM)      croak("out of memory");
    else if(status == eslENOTFOUND) croak("there is no sequence %d\n", i);
    else if(status == eslEFORMAT)   croak("error fetching sequence num %d, something wrong with SSI index?\n", i);
    else if(status != eslOK)        croak("error fetching sequence num %d\n", i);
    /* EPN, Fri Nov  6 08:44:29 2015: we allow length 0 sequences, NCBI commonly puts them in as some type 
     *                                of sequence divider for separating assembly projects, for example.
     *                                I'm leaving the following line of code here though commented out
     *                                in case it causes some unforeseen issue in the future and having it
     *                                here commented out lends some clue.
     * if(L == 0)                      croak("error fetching sequence num %d, seq length unknown\n", i); 
     */
    nres += L;
  }
  /* convert nres to a string to return (so we don't overflow an int or a long) */
  /* determine number of digits in nres */
  nres_dup = nres;
  nres_wid = 1;
  while(nres_dup >= 10) { nres_dup/=10; nres_wid++; }
  ESL_ALLOC(nres_str, sizeof(char) * (nres_wid + 1));
  snprintf(nres_str, nres_wid+1, "%" PRId64 "", nres);

  nres_str_SV = newSVpv(nres_str, nres_wid);
  free(nres_str);

  return nres_str_SV;

 ERROR: 
  croak("out of memory");
  return NULL; /* NEVER REACHED */
}
