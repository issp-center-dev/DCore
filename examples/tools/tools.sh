# Check if command CMD works
# Usage: check_command CMD
check_command () {
    which $1
    if [ $? = 1 ]; then
        echo "Command $1 not found" >&2
        exit 1
    fi
}

# Check if environment variable VAR is defined
# Usage: check_var VAR
check_var () {
    if [ -z `eval echo '$'$1` ]; then
        echo "Variable $1 not defined" >&2
        exit 1
    fi
}

# Set environment variable NUM_PROC=1 if not defined
set_num_proc () {
    echo "Setting NUM_PROC =" ${NUM_PROC:=1}
}

# Check if the previous command succeeded
check_status () {
    if [ $? = 1 ]; then
        exit 1
    fi
}
