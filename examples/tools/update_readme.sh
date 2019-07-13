#!/bin/sh

# generate "examples_info.md"
echo "### Running print_info.py"
python tools/print_info.py

# insert examples_info.md in README.md
echo ""
echo "### Updating README.md"

start="INFO TABLE START"
end="INFO TABLE END"

sed -i "/$start/a SSS" README.md  # insert 'SSS' after $start
sed -i "/$end/i EEE" README.md    # insert 'EEE' before $end
sed -i "/SSS/,/EEE/d" README.md   # delete between 'SSS' and 'EEE'
sed -i "/$start/r examples_info.md" README.md  # insert file after $start

echo " updated"
