DESCRIPTION

  This program is (for now) a simple log-odds calculator. It does not do much
  else yet, and will eventually have more utility in the prediction of U12
  introns.

INPUT

  The input for now is one of the .tsv files found in the inputs directory.
  You can specify the file you want as the first command-line argument, or omit
  it to use whichever it was compiled to use.

OUTPUT

  The output from the program is a matrix like the input showing the PWM scores
  of the given nucleotides at each position. It will also print the PWM scores
  normalized like so:
    - Values >= 0 are scaled to fall within the 50-100 range
    - Values < 0 are scaled to fall within the 0-50 range, where 0 is the
      lowest ("most negative") number
