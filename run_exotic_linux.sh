#!/usr/bin/env bash

py_commands="python3 python"
pip_commands="pip3 pip"
ver_py_min="3.6"
ver_py_max="4.0"
py_download="https://www.python.org/downloads/"
pip_download="https://bootstrap.pypa.io/get-pip.py"
pip_instructions="https://pip.pypa.io/en/stable/installing/"
py_instructions_linux="https://www.cyberciti.biz/faq/install-python-linux/"
py_runner=""
pip_runner=""

if echo "${0##*/}" | grep -q '_linux';
then
    py_download=${py_instructions_linux}
fi

# test version compatibility within predefined range
version_within_range () {
    local test_str="${1}\n${2}\n${3}"
    test_str=$(sort -V <(echo -e "${test_str}"))
    if [ "$(head -1 <(echo -e "${test_str}"))" == "${1}" ] &&
        [ "$(tail -1 <(echo -e "${test_str}"))" == "${3}" ];
    then
        return 0
    fi
    return 1
}

# test for python installation and version compatibility
for app in ${py_commands} ; do
    py_version=$(${app} --version 2>&1)
    if [ ${?} -eq 0 ];  # success
    then
        ver=$(echo -n "${py_version}" | sed 's/.*[Pp]ython \([0-9.]*\).*/\1/;' | head -1)
        version_within_range "${ver_py_min}" "${ver}" "${ver_py_max}"
        if [ $? -eq 0 ]
        then
            py_runner="${app}"
            echo "INFO: Using '${app} ${ver}'. SUCCESS!"
            break
        else
            echo "WARNING: Incompatible Python (${app}) found. ..."
        fi
    fi
done
# exit if valid python not found
if [[ -z "${py_runner}" ]];
then
    echo -e "ERROR: Incompatible or missing Python runtime. Please install\n" \
            "       Python 3.6 or above. EXITING!"
    echo -e "For more information, see ${py_download}. ...\n"
    exit 65
fi
# test for pip
for app in ${pip_commands} ; do
    pip_version=$(${app} --version 2>&1)
    if [ ${?} -eq 0 ];  # success
    then
        ver=$(echo -n "${pip_version}" | sed 's/.*[Pp]ython \([0-9.]*\).*/\1/;' | head -1)
        ver_min=$(echo -n "${ver_py_min}" | sed 's/^\([0-9]*\.[0-9]*\).*/\1/;' | head -1)
        version_within_range "${ver_min}" "${ver}" "${ver_py_max}"
        if [ $? -eq 0 ];
        then
            pip_runner="${app}"
            echo "INFO: Using '${app} ${ver}'. SUCCESS!"
            break
        else
            echo "WARNING: Incompatible Pip (${app}) found. ..."
        fi
    fi
done
# attempt install then exit on failure
if [[ -z "${pip_runner}" ]];
then
    echo "INFO: Installing 'Pip Installs Packages' (pip) for Python3. ..."
    echo "      DOWNLOADING ..."
    # install pip3
    if curl &>/dev/null ;
    then
        curl -O "${pip_download}"
    elif wget &>/dev/null ;
    then
        wget "${pip_download}"
    else
        echo -e "ERROR: Unable to download package manager. Please install\n" \
                "       Pip for Python 3. EXITING!"
        echo -e "For more information, see ${pip_instructions}. ...\n"
        exit 65
    fi
    ${py_runner} get-pip.py
    pip_runner="pip3"
    # validate installation
    if ! ${pip_runner} --version ;
    then
        echo -e "ERROR: Incompatible or missing package manager. Please install\n" \
                "       Pip for Python 3. EXITING!"
        echo -e "For more information, see ${pip_instructions}. ...\n"
        exit 65
    fi
fi
# exec commands using determinate pip
echo "INFO: Installing EXOTIC build dependencies. ..."
${pip_runner} install wheel
${pip_runner} install setuptools
echo "INFO: Installing EXOTIC core. ..."
${pip_runner} install --upgrade exotic
echo "INFO: Launching EXOTIC user interface. ..."
bash -c "exotic-gui"
