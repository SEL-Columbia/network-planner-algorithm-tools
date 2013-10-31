#!/bin/sh

# This shell script is a Makefile generator, for the quality control (QC)
# criteria we want to track automatically: it gets all list of all the
# python code files and creates a Makefile which will run them through
# pylint; finally, it generates a markdown summary (PYLINT_SCORES.md)
# of those results, timestamped with when the 'make all' command was run.

# list all the python files we want to score with pylint
pyfiles=`cd ..; ls *.py | grep -v "__init" ; ls modules/*.py | grep -v "__init"`

# Define the 'clean' rule
echo "\nclean:\n\trm *.lint *.md"

# Define the 'all' rule
echo -n "\nall:"
for pyf in $pyfiles
do
	output=`echo $pyf | tr '/' '-'`
	echo "\t$output.lint \\"
done
echo "\tPYLINT_SCORES.md\n"

# Define each individual pylint file score rule
for pyf in $pyfiles
do
	output=`echo $pyf | tr '/' '-'`
	echo "$output.lint:\n\t- pylint "../$pyf" >$output.lint 2>&1"
done

# Define the rule to create the markdown score summary file
# (use a simple grep/sed combination to pull the summary info out of
# each pylint result, though we're saving the full pylint results here
# for reference)
echo """\nPYLINT_SCORES.md:
	echo -n \"# Pylint scores as of \" > PYLINT_SCORES.md
	echo \"`date +\"%b %d %Y %H:%M:%S (UTC)\" --utc`\" >> PYLINT_SCORES.md
	echo \"<table><tr><td><b>Script</b></td><td><b>Pylint Score</b></td></tr>\" >> PYLINT_SCORES.md 
	grep \"Your code has been rated at\" *.lint | \
sed -e 's/Your code has been rated at/ /' | \
sed -e 's/modules-/modules\//' | \
sed -e 's/^/<tr><td>/' | \
sed -e 's/\.lint:/<\/td><td align=\"right\"><tt>/' | \
sed -e 's/\$\$/<\/tt><\/td><\/tr>/' >> PYLINT_SCORES.md
	echo \"</table>\" >> PYLINT_SCORES.md
"""

