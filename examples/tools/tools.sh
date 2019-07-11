check_command () {
    which $1
    if [ $? = 1 ]; then
        echo "Command $1 not found" >&2
        exit 1
    fi
}

check_var () {
    if [ -z `eval echo '$'$1` ]; then
        echo "Variable $1 not defined" >&2
        exit 1
    fi
}

set_num_proc () {
    echo "Setting NUM_PROC =" ${NUM_PROC:=1}  # set 1 if not defined
}
