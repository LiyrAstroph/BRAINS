/*
 * This code is based on the following sources, with modifications
 * 1) GNU getopt.c
 * 2) ya_getopt  - Yet another getopt (https://github.com/kubo/ya_getopt)
 *    (Copyright 2015 Kubo Takehiro <kubo@jiubao.org>)
 */

#ifndef MYGETOPT_H
#define MYGETOPT_H 1

#if defined(__cplusplus)
extern "C" {
#endif

/* For communication from `getopt' to the caller.
   When `getopt' finds an option that takes an argument,
   the argument value is returned here.
   Also, when `ordering' is RETURN_IN_ORDER,
   each non-option ARGV-element is returned here.  */

extern char *my_optarg;

/* Index in ARGV of the next element to be scanned.
  This is used for communication to and from the caller
  and for communication between successive calls to `getopt'.

  On entry to `getopt', zero means this is the first call; initialize.

  When `getopt' returns -1, this is the index of the first of the
  non-option elements that the caller should itself scan.

  Otherwise, `optind' communicates from one call to the next
  how much of ARGV has been scanned so far.  */

extern int my_optind;

/* Callers store zero here to inhibit the error message `getopt' prints
  for unrecognized options.  */

extern int my_opterr;

/* Set to an option character which was unrecognized.  */

extern int my_optopt;

/* Describe the long-named options requested by the application.
   The LONG_OPTIONS argument to getopt_long or getopt_long_only is a vector
   of `struct option' terminated by an element containing a name which is
   zero.

   The field `has_arg' is:
   no_argument		(or 0) if the option does not take an argument,
   required_argument	(or 1) if the option requires an argument,
   optional_argument 	(or 2) if the option takes an optional argument.

   If the field `flag' is not NULL, it points to a variable that is set
   to the value given in the field `val' when the option is found, but
   left unchanged if the option is not found.

   To have a long-named option do something other than set an `int' to
   a compiled-in constant, such as set a value from `optarg', set the
   option's `flag' field to zero and its `val' field to a nonzero
   value (the equivalent single-letter option character, if there is
   one).  For long options that have a zero `flag' field, `getopt'
   returns the contents of the `val' field.  */


struct my_option {
  const char *name;
  /* has_arg can't be an enum because some compilers complain about
     type mismatches in all the code that assumes it is an int.  */
  int has_arg;
  int *flag;
  int val;
};
   
/* Names for the values of the `has_arg' field of `struct option'.  */

#define	my_no_argument		    0
#define my_required_argument	1
#define my_optional_argument	2

extern int my_getopt(int argc, char * const argv[], const char *optstring);
extern int my_getopt_long(int argc, char * const argv[], const char *optstring,
                          const struct my_option *longopts, int *longindex);
extern int my_getopt_long_only(int argc, char * const argv[], const char *optstring,
                               const struct my_option *longopts, int *longindex);


#define MY_GETOPT_NO_COMPAT_MACRO
                               
#ifndef MY_GETOPT_NO_COMPAT_MACRO
#define getopt my_getopt
#define getopt_long my_getopt_long
#define getopt_long_only my_getopt_long_only
#define optarg my_optarg
#define optind my_optind
#define opterr my_opterr
#define optopt my_optopt
#define no_argument my_no_argument
#define required_argument my_required_argument
#define optional_argument my_optional_argument
#endif

#if defined(__cplusplus)
}
#endif

#endif

