#!/usr/bin/env bash

py_commands="python3 python"
pip_commands="pip3 pip"
ver_py_min="3.8"
ver_py_max="4.0"
py_download="https://www.python.org/downloads/"
pip_download="https://bootstrap.pypa.io/get-pip.py"
exotic_url="https://github.com/rzellem/EXOTIC"
test_url="https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
#test_url="https://exoplanetarchive.ipac.caltech.edu/TAP/sync?query=select%201%20from%20DUAL"
cert_instructions="https://stackoverflow.com/a/57795811/325452 or https://superuser.com/a/719047"
pip_instructions="https://pip.pypa.io/en/stable/installing/"
py_instructions_linux="https://www.cyberciti.biz/faq/install-python-linux/"
py_runner=""
pip_runner=""
test_result=""

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
    echo "ERROR: Incompatible or missing Python runtime. Please install"
    echo "       Python ${ver_py_min} or above. EXITING!"
    echo "For more information, see ${py_download}. ..."
    echo
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
    echo "INFO: Installing 'Pip Installs Packages' (pip) for Python ${ver_py_min%%.*}. ..."
    echo "      DOWNLOADING ..."
    # install pip3
    if curl --version &>/dev/null ;
    then
        curl -O "${pip_download}"
    elif wget --version &>/dev/null ;
    then
        wget "${pip_download}"
    else
        echo "ERROR: Unable to download package manager. Please install"
        echo "       Pip for Python ${ver_py_min%%.*}. EXITING!"
        echo "For more information, see ${pip_instructions}. ..."
        echo
        exit 65
    fi
    ${py_runner} get-pip.py
    pip_runner="pip3"
    # validate installation
    if ! ${pip_runner} --version ;
    then
        echo "ERROR: Incompatible or missing package manager. Please install"
        echo "       Pip for Python ${ver_py_min%%.*}. EXITING!"
        echo "For more information, see ${pip_instructions}. ..."
        echo
        exit 65
    fi
fi
echo "INFO: Validating certificate store. ..."
test_result=$(${py_runner} -u -c "import urllib.request; urllib.request.urlopen('${test_url}')" 2>&1)
if grep -q 'CERTIFICATE_VERIFY_FAILED' <<< "${test_result}" ; 
then 
    echo "ERROR: Incompatible or missing network security certificates. Please install"
    echo "       and configure the 'certifi' module from PyPi. EXITING!"
    echo "       Exa.: 'pip install --upgrade certifi' ... Then symlink, if needed."
    echo "For more information, see ${cert_instructions}. ..."
    echo
    exit 65
fi
# exec commands using determinate pip
echo "INFO: Installing EXOTIC build dependencies. ..."
${pip_runner} install setuptools
${pip_runner} install wheel
echo "INFO: Installing EXOTIC core. ..."
${pip_runner} install --upgrade exotic
echo "INFO: Launching EXOTIC user interface. ..."
if ${pip_runner} freeze | grep -iq 'exotic' ;
then
    if hash exotic-gui > /dev/null ;
    then
        bash -c "exotic-gui"
    # Ubuntu doesn't put this on path...
    else
        if [ -x ~/.local/bin/exotic-gui ] ;
        then
            bash -c "~/.local/bin/exotic-gui"
        else
            echo "ERROR: cannot find exotic-gui"
        fi
    fi
else
    echo "ERROR: Unable to launch EXOTIC, installation failed. Please verify installation"
    echo "       steps reported on screen or open a support ticket at: "
    echo "       ${exotic_url}/issues/new/choose ."
    echo "       EXITING!"
    exit 65
fi
echo "COMPLETE, EXITING!"
exit 0
