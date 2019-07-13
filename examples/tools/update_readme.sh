#!/bin/sh


# generate "examples_info.md"
python tools/print_info.py

# str=`cat examples_info.md`
# echo $str

# sed -n "/START/p" README.md

# sed -e "/START/,/END/d" README.md


# sed -e "/START/a EE" README.md

start="INFO TABLE START"
end="INFO TABLE END"

# 1: insert 'SSS' after $start and 'EEE' before $end
# 2: delete between 'SSS' and 'EEE'
# 3: insert after $start
#sed -e "/$start/a SSS" -e "/$end/i EEE" README.md \
# | sed -e "/SSS/,/EEE/d" \
# | sed -e "/$start/r examples_info.md"

# > README.md


sed -i "/$start/a SSS" README.md  # insert 'SSS' after $start
sed -i "/$end/i EEE" README.md    # insert 'EEE' before $end
sed -i "/SSS/,/EEE/d" README.md   # delete between 'SSS' and 'EEE'
sed -i "/$start/r examples_info.md" README.md  # insert file after $start


# sed -e "/START/a EE" README.md | sed -e "/START/r examples_info.md"
