
check_command () {
    which $1
    if [ $? = 1 ]; then
        echo "command $1 not found" >&2
        exit 1
    fi
}

check_var () {
    if [ -z `eval echo '$'$1` ]; then
        echo "variable $1 not defined" >&2
        exit 1
    fi
}
