#define _DEFAULT_SOURCE /* strncasecmp (Linux) */

/*
  locusDefiner.c

  Purpose
  -------
  This program implements the following bash logic as a single C program.

    1. Extract only significant SNPs
    2. Assign each SNP an interval [POS-window, POS+window]
    3. Merge overlapping intervals
    4. Keep the SNP with the smallest P value within each merged interval
       as the representative SNP

  Expected input
  --------------
  - Headered TSV.GZ
  - Header row: default columns CHR, POS, P (case-insensitive); overrides via --col-*
  - Optional: P column as -log10(P) via --p-scale neglog10
  - Sorted by chromosome then position (same order as the position column)

  Expected output
  ---------------
  Full representative sumstats row + locus_chr, locus_start, locus_end (merged interval)

  Important
  ---------
  - This is NOT physmerge-style logic that expands blocks based on the
    representative SNP
  - This is simply the union of ±window intervals around significant SNPs
*/

#include <stdio.h>      // fprintf, fopen, fclose
#include <stdlib.h>     // strtod, strtoll, exit
#include <string.h>     // strcmp, strncpy, strchr, memcpy, memcmp
#include <strings.h>      // strncasecmp
#include <stdbool.h>    // bool, true, false
#include <float.h>      // DBL_MAX
#include <stdint.h>     // int64_t
#include <math.h>       // isfinite, log10, INFINITY
#include <zlib.h>       // gzopen, gzgets, gzclose

typedef enum {
    P_SCALE_P = 0,          /* P column holds plain P in (0,1] */
    P_SCALE_NEGLOG10 = 1    /* P column holds -log10(P); convert to P internally */
} PScale;

/*
  Extremely long lines are not usually expected in ordinary summary statistics,
  so for now we set the limit to 1 MB.
  If stronger robustness is needed later, this can be replaced with a
  dynamically growing buffer.
*/
#define LINE_BUF_SIZE 1048576

/*
  Fixed-length buffer for chromosome names.
  This is enough for values such as:
  "1", "22", "X", "chr1", "chr22", "chrX", etc.
*/
#define CHR_BUF_SIZE 128

/* Upper bounds for passthrough of non-analytic sumstats columns */
#define MAX_HEADER_FIELDS 512
#define HEADER_NAME_LEN 256

/*
  Structure representing one merged locus.

  active:
    Whether this structure currently holds a valid interval

  chr:
    Chromosome name

  start, end:
    Coordinates of the current merged interval

  rps_p / rps_neglog10:
    Internal only; choose the representative SNP.

  rps_sumstat_line:
    Full representative sumstats row (no trailing newline); output appends every
    field after chr/start/end.

  rps_neglog10, rep_from_neglog10:
    When true, representative is chosen by largest -log10(P) from the column.
*/
typedef struct {
    char chr[CHR_BUF_SIZE];
    int64_t start;
    int64_t end;
    double rps_p;
    double rps_neglog10;
    bool rep_from_neglog10;
    char rps_sumstat_line[LINE_BUF_SIZE];
    bool active;
} Interval;

/*
  Print usage information
*/
static void usage(const char *prog) {
    fprintf(stderr,
        "Usage:\n"
        "  %s --input sumstats.tsv.gz --output loci.tsv --thr 5e-8 --window 500000\n"
        "      [--col-chr NAME] [--col-pos NAME] [--col-p NAME] [--p-scale p|neglog10]\n\n"
        "Header columns (first line of TSV):\n"
        "  Default: CHR, POS, P (matched case-insensitively on the full field).\n"
        "  Optional: --col-chr / --col-pos / --col-p set exact header names for each role.\n\n"
        "P column scale (--p-scale):\n"
        "  p (default): column is plain P in (0,1].\n"
        "  neglog10:    column is -log10(P). Significance uses x > -log10(thr) (no pow underflow).\n"
        "               Representative: largest x. Output: full rep row + locus_chr/start/end.\n\n"
        "Assumptions:\n"
        "  - Input is headered TSV.GZ\n"
        "  - Input is sorted by chromosome then base-pair position\n"
        "  - Logic: significant SNPs -> [POS-window, POS+window] -> merge overlaps\n",
        prog
    );
}

/*
  Small helper function to determine whether two chromosome names are equal.
  Using strcmp is fine for now.
  For further speed optimization, chromosome names could be converted
  internally to integer codes.
*/
static bool same_chr(const char *a, const char *b) {
    return strcmp(a, b) == 0;
}

/*
  Compute the left boundary.
  Whether coordinates should be treated as BED-like 0-based or 1-based
  depends on the specification, but the original bash code used:

      ($2-5e5<0?1:$2-5e5)

  so here we clamp values smaller than 1 to 1.
*/
static int64_t left_bound(int64_t pos, int64_t window) {
    int64_t s = pos - window;
    return (s < 1) ? 1 : s;
}

/*
  Initialize an Interval as "empty"
*/
static void init_interval(Interval *x) {
    x->chr[0] = '\0';
    x->start = 0;
    x->end = 0;
    x->rps_p = DBL_MAX;   // Start with a large value so smaller P values can replace it
    x->rps_neglog10 = -INFINITY;
    x->rep_from_neglog10 = false;
    x->rps_sumstat_line[0] = '\0';
    x->active = false;
}

/*
  Copy one Interval into another
*/
static void copy_interval(Interval *dst, const Interval *src) {
    strncpy(dst->chr, src->chr, CHR_BUF_SIZE - 1);
    dst->chr[CHR_BUF_SIZE - 1] = '\0';
    dst->start = src->start;
    dst->end = src->end;
    dst->rps_p = src->rps_p;
    dst->rps_neglog10 = src->rps_neglog10;
    dst->rep_from_neglog10 = src->rep_from_neglog10;
    memcpy(dst->rps_sumstat_line, src->rps_sumstat_line, sizeof(dst->rps_sumstat_line));
    dst->active = src->active;
}

/*
  Copy input line without trailing CR/LF (for storing the original sumstats row).
*/
static void copy_sumstat_line_trim(const char *src, char *dst, size_t cap) {
    if (cap == 0) {
        return;
    }
    size_t i = 0;
    while (src[i] != '\0' && src[i] != '\n' && src[i] != '\r' && i < cap - 1) {
        dst[i] = src[i];
        i++;
    }
    dst[i] = '\0';
}

/*
  Print every sumstats field from `line` (tab-separated, no trailing newline).
  No leading tab before the first field.
*/
static void fprint_all_sumstat_fields(FILE *out, const char *line) {
    const char *s = line;
    bool first = true;

    while (*s != '\0' && *s != '\n' && *s != '\r') {
        const char *fs = s;
        const char *fe = fs;
        while (*fe != '\0' && *fe != '\t' && *fe != '\n' && *fe != '\r') {
            fe++;
        }

        if (!first) {
            fputc('\t', out);
        }
        first = false;
        for (const char *q = fs; q < fe; q++) {
            fputc((unsigned char)*q, out);
        }

        if (*fe == '\t') {
            s = fe + 1;
        } else {
            break;
        }
    }
}

/*
  Print one interval to the output file
*/
static void print_interval(FILE *out, const Interval *x) {
    fprint_all_sumstat_fields(out, x->rps_sumstat_line);
    fprintf(out, "\t%s\t%lld\t%lld\n",
            x->chr,
            (long long)x->start,
            (long long)x->end);
}

/*
  True if the header field [fs, fs+len) equals `name` (case-insensitive, full field).
*/
static bool field_name_match(const char *fs, size_t len, const char *name) {
    size_t nlen = strlen(name);
    return len > 0 && len == nlen && strncasecmp(fs, name, len) == 0;
}

/*
  Resolve column indices from the header line.

  col_chr / col_pos / col_p:
    If NULL, the default name is CHR, POS, or P respectively.
    Otherwise the given string is the exact header label to match (case-insensitive).

  Duplicate labels in the file: the last matching column wins for that role.

  field_names:
    Filled with each header cell (column names), up to MAX_HEADER_FIELDS.

  n_fields_out:
    Total number of columns in the header row.

  Return value:
    0  : success
   -1  : one or more required columns are missing, or too many columns
*/
static int parse_header(const char *line, int *chr_idx, int *pos_idx, int *p_idx,
                        const char *col_chr, const char *col_pos, const char *col_p,
                        char (*field_names)[HEADER_NAME_LEN], int *n_fields_out) {
    const char *name_chr = col_chr ? col_chr : "CHR";
    const char *name_pos = col_pos ? col_pos : "POS";
    const char *name_p = col_p ? col_p : "P";

    *chr_idx = *pos_idx = *p_idx = -1;
    int idx = 0;
    const char *s = line;

    while (*s != '\0' && *s != '\n' && *s != '\r') {
        if (idx >= MAX_HEADER_FIELDS) {
            return -1;
        }

        const char *fs = s;
        const char *fe = fs;
        while (*fe != '\0' && *fe != '\t' && *fe != '\n' && *fe != '\r') {
            fe++;
        }
        size_t len = (size_t)(fe - fs);

        const char *nfs = fs;
        size_t nlen = len;
        if (idx == 0 && len >= 3 &&
            (unsigned char)fs[0] == 0xEF &&
            (unsigned char)fs[1] == 0xBB &&
            (unsigned char)fs[2] == 0xBF) {
            nfs += 3;
            nlen -= 3;
        }

        {
            size_t copy = nlen;
            if (copy >= HEADER_NAME_LEN) {
                copy = HEADER_NAME_LEN - 1;
            }
            memcpy(field_names[idx], nfs, copy);
            field_names[idx][copy] = '\0';
        }

        if (field_name_match(nfs, nlen, name_chr)) {
            *chr_idx = idx;
        }
        if (field_name_match(nfs, nlen, name_pos)) {
            *pos_idx = idx;
        }
        if (field_name_match(nfs, nlen, name_p)) {
            *p_idx = idx;
        }

        if (*fe == '\t') {
            s = fe + 1;
            idx++;
        } else {
            break;
        }
    }

    *n_fields_out = idx + 1;
    return (*chr_idx >= 0 && *pos_idx >= 0 && *p_idx >= 0) ? 0 : -1;
}

/*
  Extract CHR / POS / P by scanning tabs once. Stops after the last required
  column once all three fields are filled. Temporarily NUL-terminates fields
  for strtoll/strtod and restores the delimiter before continuing.

  No per-line copy into a second LINE_BUF_SIZE buffer: the old strncpy(work, line,
  LINE_BUF_SIZE-1) pattern pads to nearly 1 MiB on many libcs for short rows,
  which dominated runtime on large files.

  Return value:
    0  : success
   -1  : parse failure
*/
static int parse_record(char *line, int chr_idx, int pos_idx, int p_idx,
                        char *chr_out, int64_t *pos_out, double *p_out) {
    bool got_chr = false;
    bool got_pos = false;
    bool got_p = false;
    int max_i = chr_idx;
    if (pos_idx > max_i) {
        max_i = pos_idx;
    }
    if (p_idx > max_i) {
        max_i = p_idx;
    }

    char *s = line;
    int idx = 0;

    for (;;) {
        char *fs = s;
        char *fe = fs;
        while (*fe != '\0' && *fe != '\t' && *fe != '\n' && *fe != '\r') {
            fe++;
        }

        if (idx == chr_idx) {
            size_t len = (size_t)(fe - fs);
            if (len >= CHR_BUF_SIZE) {
                len = CHR_BUF_SIZE - 1;
            }
            memcpy(chr_out, fs, len);
            chr_out[len] = '\0';
            got_chr = true;
        } else if (idx == pos_idx) {
            char saved = *fe;
            *fe = '\0';
            char *endptr = NULL;
            long long v = strtoll(fs, &endptr, 10);
            if (endptr == fs || *endptr != '\0') {
                *fe = saved;
                return -1;
            }
            *fe = saved;
            *pos_out = (int64_t)v;
            got_pos = true;
        } else if (idx == p_idx) {
            char saved = *fe;
            *fe = '\0';
            char *endptr = NULL;
            double v = strtod(fs, &endptr);
            if (endptr == fs || *endptr != '\0') {
                *fe = saved;
                return -1;
            }
            *fe = saved;
            *p_out = v;
            got_p = true;
        }

        if (idx >= max_i && got_chr && got_pos && got_p) {
            return 0;
        }

        if (*fe == '\t') {
            s = fe + 1;
            idx++;
        } else {
            break;
        }
    }

    return (got_chr && got_pos && got_p) ? 0 : -1;
}

/*
  Feed one significant SNP into the current merged interval.

  next:
    A new interval created from this SNP alone:
      [pos-window, pos+window]

  cur:
    The currently accumulated merged interval

  Cases:
    1. If cur is empty, start with next
    2. If the chromosome changes, or next does not overlap cur,
       output cur and start a new interval with next
    3. If they overlap, merge them and update the representative SNP if needed
*/
static void consume_sig_snp(FILE *out, Interval *cur,
                            const char *chr, int64_t pos, double p_lin,
                            int64_t window, bool rep_from_neglog10, double neglog10,
                            const char *sumstat_line) {
    Interval next;

    strncpy(next.chr, chr, CHR_BUF_SIZE - 1);
    next.chr[CHR_BUF_SIZE - 1] = '\0';
    next.start = left_bound(pos, window);
    next.end = pos + window;
    next.rps_p = p_lin;
    next.rps_neglog10 = neglog10;
    next.rep_from_neglog10 = rep_from_neglog10;
    copy_sumstat_line_trim(sumstat_line, next.rps_sumstat_line, sizeof(next.rps_sumstat_line));
    next.active = true;

    /* If nothing has been accumulated yet, start with next */
    if (!cur->active) {
        copy_interval(cur, &next);
        return;
    }

    /*
      If the chromosome changed, or next does not overlap the current interval:
      the current interval is finalized, so output it and start a new one.

      Non-overlap is defined here as:
          next.start > cur->end

      Therefore, if
          next.start <= cur->end
      then the intervals are merged.
    */
    if (!same_chr(cur->chr, next.chr) || next.start > cur->end) {
        print_interval(out, cur);
        copy_interval(cur, &next);
        return;
    }

    /*
      Reaching here means cur and next overlap on the same chromosome.
      Therefore, take their union.
    */
    if (next.end > cur->end) {
        cur->end = next.end;
    }
    if (next.start < cur->start) {
        cur->start = next.start;
    }

    /*
      Representative SNP: smallest linear P, i.e. largest -log10(P) when using
      neglog10-scale input (avoids comparing underflowed pow values).
    */
    if (rep_from_neglog10) {
        if (neglog10 > cur->rps_neglog10) {
            cur->rps_neglog10 = neglog10;
            copy_sumstat_line_trim(sumstat_line, cur->rps_sumstat_line,
                                   sizeof(cur->rps_sumstat_line));
        }
    } else {
        if (p_lin < cur->rps_p) {
            cur->rps_p = p_lin;
            copy_sumstat_line_trim(sumstat_line, cur->rps_sumstat_line,
                                   sizeof(cur->rps_sumstat_line));
        }
    }
}

static bool col_opt_invalid(const char *s) {
    return s == NULL || s[0] == '\0' ||
           strchr(s, '\t') != NULL || strchr(s, '\n') != NULL || strchr(s, '\r') != NULL;
}

int main(int argc, char **argv) {
    const char *input_path = NULL;
    const char *output_path = NULL;
    double thr = 5e-8;
    int64_t window = 500000;
    const char *col_chr = NULL;
    const char *col_pos = NULL;
    const char *col_p = NULL;
    PScale p_scale = P_SCALE_P;

    /*
      Read command-line arguments.
      This is kept simple and does not use getopt.
    */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--input") == 0 && i + 1 < argc) {
            input_path = argv[++i];
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output_path = argv[++i];
        } else if (strcmp(argv[i], "--thr") == 0 && i + 1 < argc) {
            thr = strtod(argv[++i], NULL);
        } else if (strcmp(argv[i], "--window") == 0 && i + 1 < argc) {
            window = (int64_t)strtoll(argv[++i], NULL, 10);
        } else if (strcmp(argv[i], "--col-chr") == 0 && i + 1 < argc) {
            col_chr = argv[++i];
        } else if (strcmp(argv[i], "--col-pos") == 0 && i + 1 < argc) {
            col_pos = argv[++i];
        } else if (strcmp(argv[i], "--col-p") == 0 && i + 1 < argc) {
            col_p = argv[++i];
        } else if (strcmp(argv[i], "--p-scale") == 0 && i + 1 < argc) {
            const char *v = argv[++i];
            if (strcmp(v, "p") == 0) {
                p_scale = P_SCALE_P;
            } else if (strcmp(v, "neglog10") == 0) {
                p_scale = P_SCALE_NEGLOG10;
            } else {
                fprintf(stderr, "ERROR: --p-scale must be `p` or `neglog10`.\n");
                return 1;
            }
        } else {
            usage(argv[0]);
            return 1;
        }
    }

    if (input_path == NULL || output_path == NULL || thr <= 0.0 || window <= 0) {
        usage(argv[0]);
        return 1;
    }

    if (col_opt_invalid(col_chr) && col_chr != NULL) {
        fprintf(stderr, "ERROR: --col-chr requires a non-empty name without tab/newline.\n");
        return 1;
    }
    if (col_opt_invalid(col_pos) && col_pos != NULL) {
        fprintf(stderr, "ERROR: --col-pos requires a non-empty name without tab/newline.\n");
        return 1;
    }
    if (col_opt_invalid(col_p) && col_p != NULL) {
        fprintf(stderr, "ERROR: --col-p requires a non-empty name without tab/newline.\n");
        return 1;
    }

    /*
      Open gzip-compressed input
    */
    gzFile gz = gzopen(input_path, "rb");
    if (!gz) {
        fprintf(stderr, "ERROR: failed to open input: %s\n", input_path);
        return 1;
    }

    /*
      Open output file
    */
    FILE *out = fopen(output_path, "w");
    if (!out) {
        fprintf(stderr, "ERROR: failed to open output: %s\n", output_path);
        gzclose(gz);
        return 1;
    }

    /*
      Use a larger output buffer to slightly reduce the overhead
      of repeated fprintf calls.
    */
    setvbuf(out, NULL, _IOFBF, 1 << 20);

    char line[LINE_BUF_SIZE];
    char field_names[MAX_HEADER_FIELDS][HEADER_NAME_LEN];

    /*
      Read the first line = header
    */
    if (gzgets(gz, line, LINE_BUF_SIZE) == NULL) {
        fprintf(stderr, "ERROR: empty input.\n");
        gzclose(gz);
        fclose(out);
        return 1;
    }

    int chr_idx, pos_idx, p_idx;
    int n_fields = 0;
    if (parse_header(line, &chr_idx, &pos_idx, &p_idx, col_chr, col_pos, col_p,
                      field_names, &n_fields) != 0) {
        fprintf(stderr,
                "ERROR: header must include chromosome, position, and P-value columns "
                "(default names CHR, POS, P; case-insensitive), or names set by "
                "--col-chr / --col-pos / --col-p; at most %d columns.\n",
                MAX_HEADER_FIELDS);
        gzclose(gz);
        fclose(out);
        return 1;
    }

    /*
      Output header: full input header + merged-locus suffix
    */
    for (int hi = 0; hi < n_fields; hi++) {
        if (hi > 0) {
            fputc('\t', out);
        }
        fputs(field_names[hi], out);
    }
    fprintf(out, "\tlocus_chr\tlocus_start\tlocus_end\n");

    /*
      Current merged interval being accumulated
    */
    Interval cur;
    init_interval(&cur);

    /*
      Variables for checking input sort order
      - Reset last_pos when chromosome changes
      - Raise an error if positions go backwards within a chromosome
    */
    char last_chr[CHR_BUF_SIZE] = "";
    bool seen_any = false;
    int64_t last_pos = -1;

    double thr_neglog10 = 0.0;
    if (p_scale == P_SCALE_NEGLOG10) {
        thr_neglog10 = -log10(thr);
        if (!isfinite(thr_neglog10)) {
            fprintf(stderr, "ERROR: --thr is invalid for --p-scale neglog10.\n");
            gzclose(gz);
            fclose(out);
            return 1;
        }
    }

    /*
      Main loop:
      read one line at a time and process it immediately
    */
    while (gzgets(gz, line, LINE_BUF_SIZE) != NULL) {
        /*
          Simple detection of overlong lines.
          If no newline is present and we are not at EOF,
          assume the line exceeds the buffer.
        */
        if (strchr(line, '\n') == NULL && !gzeof(gz)) {
            fprintf(stderr, "ERROR: line exceeds buffer size (%d).\n", LINE_BUF_SIZE);
            gzclose(gz);
            fclose(out);
            return 1;
        }

        char chr[CHR_BUF_SIZE];
        int64_t pos = 0;
        double p_col = 1.0;

        /*
          Skip malformed lines.
          In production, this could be made stricter and turned into an error.
        */
        if (parse_record(line, chr_idx, pos_idx, p_idx, chr, &pos, &p_col) != 0) {
            continue;
        }

        if (p_scale == P_SCALE_NEGLOG10) {
            if (!isfinite(p_col) || p_col < 0.0) {
                continue;
            }
        } else {
            if (!(p_col > 0.0 && p_col <= 1.0)) {
                continue;
            }
        }

        /*
          Check sort order.
          Since this program performs streaming merging,
          sorted input is essential.
        */
        if (!seen_any) {
            strncpy(last_chr, chr, CHR_BUF_SIZE - 1);
            last_chr[CHR_BUF_SIZE - 1] = '\0';
            last_pos = pos;
            seen_any = true;
        } else {
            if (!same_chr(chr, last_chr)) {
                /*
                  If the chromosome changes, finalize and output the current
                  interval from the previous chromosome, then reset.
                */
                if (cur.active) {
                    print_interval(out, &cur);
                    init_interval(&cur);
                }

                strncpy(last_chr, chr, CHR_BUF_SIZE - 1);
                last_chr[CHR_BUF_SIZE - 1] = '\0';
                last_pos = pos;
            } else {
                /*
                  Raise an error if positions go backwards within the same chromosome.
                */
                if (pos < last_pos) {
                    fprintf(stderr,
                            "ERROR: input is not sorted within chromosome %s (%lld after %lld).\n",
                            chr, (long long)pos, (long long)last_pos);
                    gzclose(gz);
                    fclose(out);
                    return 1;
                }
                last_pos = pos;
            }
        }

        /*
          Only significant SNPs are passed into the merge logic.
          neglog10: compare x to -log10(thr) so pow(10,-x) never underflows here.
        */
        if (p_scale == P_SCALE_NEGLOG10) {
            if (p_col > thr_neglog10) {
                consume_sig_snp(out, &cur, chr, pos, 0.0, window, true, p_col, line);
            }
        } else {
            if (p_col < thr) {
                consume_sig_snp(out, &cur, chr, pos, p_col, window, false, 0.0, line);
            }
        }
    }

    /*
      After the loop, output the final remaining interval if present.
    */
    if (cur.active) {
        print_interval(out, &cur);
    }

    gzclose(gz);
    fclose(out);
    return 0;
}
