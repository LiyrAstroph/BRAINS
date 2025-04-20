/*
 * This code is based on the following sources, with modifications
 * 1) GNU getopt.c
 * 2) ya_getopt  - Yet another getopt (https://github.com/kubo/ya_getopt)
 *    (Copyright 2015 Kubo Takehiro <kubo@jiubao.org>)
 */

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mygetopt.h"

/* For communication from `getopt' to the caller.
   When `getopt' finds an option that takes an argument,
   the argument value is returned here.
   Also, when `ordering' is RETURN_IN_ORDER,
   each non-option ARGV-element is returned here.  */

char *my_optarg = NULL;

/* Index in ARGV of the next element to be scanned.
  This is used for communication to and from the caller
  and for communication between successive calls to `getopt'.

  On entry to `getopt', zero means this is the first call; initialize.

  When `getopt' returns -1, this is the index of the first of the
  non-option elements that the caller should itself scan.

  Otherwise, `optind' communicates from one call to the next
  how much of ARGV has been scanned so far.  */

/* 1003.2 says this must be 1 before any call.  */
int my_optind = 1;

/* Formerly, initialization of getopt depended on optind==0, which
  causes problems with re-calling getopt as programs generally don't
  know that. */

int __getopt_initialized = 0;

/* The next char to be scanned in the option-element
   in which the last option character we returned was found.
   This allows us to pick up the scan where we left off.

   If this is zero, or a null string, it means resume the scan
   by advancing to the next ARGV-element.  */

static char *my_optnext = NULL;

/* Callers store zero here to inhibit the error message
   for unrecognized options.  */

int my_opterr = 1;

/* Set to an option character which was unrecognized.
   This must be initialized on some systems to avoid linking in the
   system's own getopt implementation.  */

int my_optopt = '?';

/* posixly_correct = 1 means not permuting argv 
 */
static int posixly_correct = -1;

/* whether handle non-option arguments
 */
static int handle_nonopt_argv = 0;

static enum
{
  MY_REQUIRE_ORDER, MY_PERMUTE, MY_RETURN_IN_ORDER
} my_ordering;

/* Describe the part of ARGV that contains non-options that have
   been skipped.  `first_nonopt' is the index in ARGV of the first of them;
   `last_nonopt' is the index after the last of them.  
   Valid only when PERMUTE */

static int my_first_nonopt;
static int my_last_nonopt;

/* Describe how to deal with options that follow non-option ARGV-elements.

  If the caller did not specify anything,
  the default is REQUIRE_ORDER if the environment variable
  POSIXLY_CORRECT is defined, PERMUTE otherwise.

  REQUIRE_ORDER means don't recognize them as options;
  stop option processing when the first non-option is seen.
  This is what Unix does.
  This mode of operation is selected by either setting the environment
  variable POSIXLY_CORRECT, or using `+' as the first character
  of the list of option characters.

  PERMUTE is the default.  We permute the contents of ARGV as we scan,
  so that eventually all the non-options are at the end.  This allows options
  to be given in any order, even with programs that were not written to
  expect this.

  RETURN_IN_ORDER is an option available to programs that were written
  to expect options and other ARGV-elements in any order and that care about
  the ordering of the two.  We describe each non-option ARGV-element
  as if it were the argument of an option with character code 1.
  Using `-' as the first character of the list of option characters
  selects this mode of operation.

  The special argument `--' forces an end of option-scanning regardless
  of the value of `ordering'.  In the case of RETURN_IN_ORDER, only
  `--' can cause `getopt' to return -1 with `optind' != ARGC.  */

static void my_getopt_error(const char *optstring, const char *format, ...);
static void check_gnu_extension(const char *optstring);
static int my_getopt_internal(int argc, char * const argv[], const char *optstring, 
                              const struct my_option *longopts, int *longindex, int long_only);
static int my_getopt_shortopts(int argc, char * const argv[], const char *optstring, int long_only);
static int my_getopt_longopts(int argc, char * const argv[], char *arg, const char *optstring, 
                              const struct my_option *longopts, int *longindex, int *long_only_flag);

/* print error to stderr */
static void my_getopt_error(const char *optstring, const char *format, ...)
{
  if (my_opterr && optstring[0] != ':') 
  {
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
  }
}

/* check gnu extension */
static void check_gnu_extension(const char *optstring)
{
  /* "+" means not permuting argv */
  if (optstring[0] == '+' || getenv("POSIXLY_CORRECT") != NULL) 
  {
    posixly_correct = 1;
    my_ordering = MY_REQUIRE_ORDER;
  } 
  else 
  {
    posixly_correct = 0;
    my_ordering = MY_PERMUTE;
  }

  /* "-" means keeping argv in order as input */
  if (optstring[0] == '-') 
  {
    handle_nonopt_argv = 1;
    my_ordering = MY_RETURN_IN_ORDER;
  } 
  else 
  {
    handle_nonopt_argv = 0;
    my_ordering = MY_PERMUTE;
  }
}

/* check whether arg is an option like "-x" 
 * a single "-" character is not an option
 */
static int is_option(const char *arg)
{
  return arg[0] == '-' && arg[1] != '\0';
}

/* handle only short options */
int my_getopt(int argc, char * const argv[], const char *optstring)
{
  return my_getopt_internal(argc, argv, optstring, NULL, NULL, 0);
}

/* handle both short and long options */
int my_getopt_long(int argc, char * const argv[], const char *optstring, 
                   const struct my_option *longopts, int *longindex)
{
  return my_getopt_internal(argc, argv, optstring, longopts, longindex, 0);
}

/* handle only long options */
int my_getopt_long_only(int argc, char * const argv[], const char *optstring, 
                        const struct my_option *longopts, int *longindex)
{
  return my_getopt_internal(argc, argv, optstring, longopts, longindex, 1);
}

/* Scan elements of ARGV (whose length is ARGC) for option characters
   given in OPTSTRING.

   If an element of ARGV starts with '-', and is not exactly "-" or "--",
   then it is an option element.  The characters of this element
   (aside from the initial '-') are option characters.  If `getopt'
   is called repeatedly, it returns successively each of the option characters
   from each of the option elements.

   If `getopt' finds another option character, it returns that character,
   updating `optind' and `nextchar' so that the next call to `getopt' can
   resume the scan with the following option character or ARGV-element.

   If there are no more option characters, `getopt' returns -1.
   Then `optind' is the index in ARGV of the first ARGV-element
   that is not an option.  (The ARGV-elements have been permuted
   so that those that are not options now come last.)

   OPTSTRING is a string containing the legitimate option characters.
   If an option character is seen that is not listed in OPTSTRING,
   return '?' after printing an error message.  If you set `opterr' to
   zero, the error message is suppressed but we still return '?'.

   If a char in OPTSTRING is followed by a colon, that means it wants an arg,
   so the following text in the same ARGV-element, or the text of the following
   ARGV-element, is returned in `optarg'.  Two colons mean an option that
   wants an optional arg; if there is text in the current ARGV-element,
   it is returned in `optarg', otherwise `optarg' is set to zero.

   If OPTSTRING starts with `-' or `+', it requests different methods of
   handling the non-option ARGV-elements.
   See the comments about RETURN_IN_ORDER and REQUIRE_ORDER, above.

   Long-named options begin with `--' instead of `-'.
   Their names may be abbreviated as long as the abbreviation is unique
   or is an exact match for some defined option.  If they have an
   argument, it follows the option name in the same ARGV-element, separated
   from the option name by a `=', or else the in next ARGV-element.
   When `getopt' finds a long-named option, it returns 0 if that option's
   `flag' field is nonzero, the value of the option's `val' field
   if the `flag' field is zero.

   The elements of ARGV aren't really const, because we permute them.
   But we pretend they're const in the prototype to be compatible
   with other systems.

   LONGOPTS is a vector of `struct option' terminated by an
   element containing a name which is zero.

   LONGIND returns the index in LONGOPT of the long-named option found.
   It is only valid when a long-named option has been found by the most
   recent call.

   If LONG_ONLY is nonzero, '-' as well as '--' can introduce
   long-named options.  */
static int my_getopt_internal(int argc, char * const argv[], const char *optstring, 
                              const struct my_option *longopts, int *longindex, int long_only)
{
  /* static start and end are initiated to 0 
   * start points to the current non-option argument
   * end points to the next option argument
   * 
   * start = 0 means no need to permute
   */
  static int start, end;

  if (my_optopt == '?') 
  {
    my_optopt = 0;
  }

  if (posixly_correct == -1) 
  {
    check_gnu_extension(optstring);
  }
  
  /* reset getopt, only check gnu extension when reset */
  if (my_optind == 0) 
  {
    check_gnu_extension(optstring);
    my_optind = 1;
    my_optnext = NULL;
    start = end = 0;
    my_first_nonopt = my_last_nonopt = 1; /* 0th argument is always program name */
  }

  /* skip the beginning "+" or "-" */
  switch (optstring[0]) 
  {
    case '+':
    case '-':
      optstring++;
  }
  
  /* advace scan and permute argv */
  if (my_optnext == NULL && start != 0) 
  {
    int last_pos = my_optind - 1; 
   
    /* because non-option arguments are placed at the end, update optind accordingly */
    my_optind -= end - start;
    if (my_optind <= 0) 
    {
      my_optind = 1;
    }

    /* place the non-option arguments at the end 
     * when end == start, non-option arguments are moved to the end, stop loop.
     */
    my_last_nonopt = end;
    while (start < end--) 
    {
      int i;
      char *arg = argv[end];

      for (i = end; i < last_pos; i++) 
      {
        ((char **)argv)[i] = argv[i + 1];
      }
      ((char const **)argv)[i] = arg;
      last_pos--;

      /* store the index of the first non-option argument */
      my_first_nonopt = i;
    }
    start = 0; 
  }

  if (my_optind >= argc) 
  {
    my_optarg = NULL;
    return -1;
  }

  if (my_optnext == NULL) 
  {
    const char *arg = argv[my_optind];
    /* arg is not an option */
    if (!is_option(arg)) 
    {
      if (handle_nonopt_argv) /* keep argv in order, 
                               * non-option arguments are handled as those of an option "1"*/
      {
        my_optarg = argv[my_optind++];
        start = 0;
        return 1;
      } 
      else if (posixly_correct) /* non-option argument, stop to parse */
      {
        my_optarg = NULL;
        return -1;
      } 
      else /* search for next option argument */
      {
        int i;

        start = my_optind;
        for (i = my_optind + 1; i < argc; i++) 
        {
          if (is_option(argv[i])) 
          {
            end = i;
            break;
          }
        }
        /* no option arguments left, return */
        if (i == argc) 
        {
          my_optarg = NULL;
          return -1;
        }
        /* point to the found option argument */
        my_optind = i;
        arg = argv[my_optind];
      }
    }
    
    /* optstring is forced to end */
    if (strcmp(arg, "--") == 0) 
    {
      my_optind++;
      return -1;
    }
    /* arg is a long option */
    if (longopts != NULL && arg[1] == '-') 
    {
      return my_getopt_longopts(argc, argv, argv[my_optind] + 2, optstring, longopts, longindex, NULL);
    }
  }

  if (my_optnext == NULL) 
  {
    my_optnext = argv[my_optind] + 1;
  }

  if (long_only) 
  {
    int long_only_flag = 0;
    int rv = my_getopt_longopts(argc, argv, my_optnext, optstring, longopts, longindex, &long_only_flag);
    if (!long_only_flag) 
    {
      my_optnext = NULL;
      return rv;
    }
  }

  return my_getopt_shortopts(argc, argv, optstring, long_only);
}

/* handle short opts */
static int my_getopt_shortopts(int argc, char * const argv[], const char *optstring, int long_only)
{
  int opt = *my_optnext;
  /* location of the opt */
  const char *os = strchr(optstring, opt);
  
  /* opt is not in optstring */
  if (os == NULL)
  {
    my_optarg = NULL;
    if (long_only) 
    {
      my_getopt_error(optstring, "%s: unrecognized option '-%s'\n", argv[0], my_optnext);
      my_optind++;
      my_optnext = NULL;
    } 
    else 
    {
      my_optopt = opt;
      my_getopt_error(optstring, "%s: invalid option -- '%c'\n", argv[0], opt);

      /* ascii number of '\0' character is 0 */
      if (*(++my_optnext) == 0) 
      {
        my_optind++;
        my_optnext = NULL;
      }
    }
    return '?';
  }
  /* opt has an argument */
  if (os[1] == ':') 
  {
    /* no argument follows: 
     * if optional argument, OK.
     * if not, return error */
    if (my_optnext[1] == 0) 
    {
      my_optind++;
      my_optnext = NULL;
      if (os[2] == ':') 
      {
        /* optional argument */
        my_optarg = NULL;
      } 
      else 
      {
        if (my_optind == argc) 
        {
          my_optarg = NULL;
          my_optopt = opt;
          my_getopt_error(optstring, "%s: option requires an argument -- '%c'\n", argv[0], opt);
          if (optstring[0] == ':') 
          {
            return ':';
          } 
          else 
          {
            return '?';
          }
        }
        my_optarg = argv[my_optind];
        my_optind++;
      }
    } 
    else 
    {
      my_optarg = my_optnext + 1;
      my_optind++;
    }
    my_optnext = NULL;
  } 
  else /* if not have an argument */
  {
    my_optarg = NULL;
    if (my_optnext[1] == 0) 
    {
      my_optnext = NULL;
      my_optind++;
    } 
    else 
    {
      my_optnext++;
    }
  }
  return opt;
}

/* handle long opts */
static int my_getopt_longopts(int argc, char * const argv[], char *arg, const char *optstring, 
                              const struct my_option *longopts, int *longindex, int *long_only_flag)
{
  char *val = NULL;
  const struct my_option *opt;
  size_t namelen;
  int idx;

  /* loop over all options */
  for (idx = 0; longopts[idx].name != NULL; idx++) 
  {
    opt = &longopts[idx];
    namelen = strlen(opt->name);
    if (strncmp(arg, opt->name, namelen) == 0) 
    {
      switch (arg[namelen]) 
      {
        case '\0':
          switch (opt->has_arg) 
          {
            case my_required_argument:
              my_optind++;
              if (my_optind == argc) 
              {
                my_optarg = NULL;
                my_optopt = opt->val;
                my_getopt_error(optstring, "%s: option '--%s' requires an argument\n", argv[0], opt->name);
                if (optstring[0] == ':') 
                {
                    return ':';
                } else 
                {
                    return '?';
                }
              }
              val = argv[my_optind];
              break;
          }
          goto found;

        case '=':
          if (opt->has_arg == my_no_argument) 
          {
            const char *hyphens = (argv[my_optind][1] == '-') ? "--" : "-";

            my_optind++;
            my_optarg = NULL;
            my_optopt = opt->val;
            my_getopt_error(optstring, "%s: option '%s%s' doesn't allow an argument\n", argv[0], hyphens, opt->name);
            return '?';
          }
          val = arg + namelen + 1;
          goto found;
      }
    }
  }
  if (long_only_flag) 
  {
    *long_only_flag = 1;
  } 
  else 
  {
    my_getopt_error(optstring, "%s: unrecognized option '%s'\n", argv[0], argv[my_optind]);
    my_optind++;
  }
  return '?';

found:
  my_optarg = val;
  my_optind++;
  if (opt->flag) 
  {
    *opt->flag = opt->val;
  }
  if (longindex) 
  {
    *longindex = idx;
  }
  return opt->flag ? 0 : opt->val;
}

/*===========================================================*/
#ifdef MYGETOPT_TEST

static void print_opts(int opt, int argc, char **argv)
{
  int i;

  if (isprint(opt)) {
      fprintf(stderr, "retval = '%c', ", opt);
  } else {
      fprintf(stderr, "retval = %3d, ", opt);
  }
  fprintf(stderr, "optind = %d, ", my_optind);
  fprintf(stderr, "opterr = %d, ", my_opterr);
  if (isprint(my_optopt)) {
      fprintf(stderr, "optopt = '%c', ", my_optopt);
  } else {
      fprintf(stderr, "optopt = %3d, ", my_optopt);
  }
  if (my_optarg == NULL) {
      fprintf(stderr, "optarg = (null)\n");
  } else {
      fprintf(stderr, "optarg = \"%s\"\n", my_optarg);
  }
  for (i = 0; i < argc; i++) {
      fprintf(stderr, "'%s'%c", argv[i], (i + 1 == argc) ? '\n' : ' ');
  }
}

void test_short(int argc, char **argv)
{
  int c;
  int digit_optind = 0;

  while (1)
  {
    int this_option_optind = my_optind ? my_optind : 1;

    c = my_getopt (argc, argv, "abc:d:0123456789");
    print_opts(c, argc, argv);

    if (c == -1)
      break;
    
    switch (c)
    {
      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        if (digit_optind != 0 && digit_optind != this_option_optind)
          printf ("digits occur in two different argv-elements.\n");
        digit_optind = this_option_optind;
        printf ("option %c\n", c);
        break;

      case 'a':
        printf ("option a\n");
        break;

      case 'b':
        printf ("option b\n");
        break;

      case 'c':
        printf ("option c with value `%s'\n", my_optarg);
        break;

      case '?':
        break;
      
      case 1:
        printf("non-option argument '%s'\n", my_optarg);
        break;

      default:
        printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }

  if (my_optind < argc)
  {
    printf ("non-option ARGV-elements: ");
    while (my_optind < argc)
      printf ("%s ", argv[my_optind++]);
    printf ("\n");
  }
  
  /* in the end, argv will be reordered. The non-option args are placed after the option args */
  int i;
  for(i=0; i<argc; i++)
  {
    printf("%s\n", argv[i]);
  }

}

void test_long(int argc, char **argv)
{
  int c;
  int digit_optind = 0;

  while (1)
  {
    int this_option_optind = my_optind ? my_optind : 1;
    int option_index = 0;
    static struct my_option long_options[] =
    {
      {"add", 1, 0, 0},
      {"append", 0, 0, 0},
      {"delete", 1, 0, 0},
      {"verbose", 0, 0, 0},
      {"create", 0, 0, 0},
      {"file", 1, 0, 0},
      {0, 0, 0, 0}
    };

    c = my_getopt_long (argc, argv, "abc:d:0123456789", long_options, &option_index);

    if (c == -1)
      break;
    
    switch (c)
    {
      case 0:
        printf ("option %s", long_options[option_index].name);
        if (my_optarg)
          printf (" with arg %s", my_optarg);
        printf ("\n");
        break;

      case '0':
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
      case '9':
        if (digit_optind != 0 && digit_optind != this_option_optind)
          printf ("digits occur in two different argv-elements.\n");
        digit_optind = this_option_optind;
        printf ("option %c\n", c);
        break;

      case 'a':
        printf ("option a\n");
        break;

      case 'b':
        printf ("option b\n");
        break;

      case 'c':
        printf ("option c with value '%s'\n", my_optarg);
        break;

      case 'd':
        printf ("option d with value '%s'\n", my_optarg);
        break;

      case '?':
        break;
      
      case 1:
        printf("non-option argument '%s'\n", my_optarg);
        break;

      default:
        printf ("?? getopt returned character code 0%o ??\n", c);
    }
  }

  if (my_optind < argc)
  {
    printf ("non-option ARGV-elements: ");
    while (my_optind < argc)
      printf ("%s ", argv[my_optind++]);
    printf ("\n");
  }

  /* in the end, argv will be reordered. The non-option args are placed after the option args */
  int i;
  for(i=0; i<argc; i++)
  {
    printf("%s\n", argv[i]);
  }
}

char **copy_argv(int argc, char **argv)
{
  int i;
  char **argv_copy;
  argv_copy = (char **)malloc(argc*sizeof(char *));
  for(i=0; i<argc; i++)
  {
    argv_copy[i] = (char *)malloc((strlen(argv[i])+1)*sizeof(char));
    strcpy(argv_copy[i], argv[i]);
  }
  return argv_copy;
}

void free_argv(int argc, char **argv)
{
  int i;
  for(i=0; i<argc; i++)
  {
    free(argv[i]);
  }
  free(argv);
}

int main(int argc, char **argv)
{
  char **argv_copy;
  printf("=============================\n");
  argv_copy = copy_argv(argc, argv);
  test_short(argc, argv_copy);
  free_argv(argc, argv_copy);
  printf("=============================\n");
  my_optind = 0;
  argv_copy = copy_argv(argc, argv);
  test_long(argc, argv);
  free_argv(argc, argv_copy);

  printf("first-nonopt: %d %d\n", my_first_nonopt, my_last_nonopt);
}

#endif
/*===========================================================*/
