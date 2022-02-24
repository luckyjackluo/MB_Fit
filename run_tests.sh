#!/bin/bash

export MBFIT_HOME=$PWD


if [ $# -eq 0 ]
then
    time=$(date '+%Y-%m-%d-%H-%M-%S')
    mkdir -p test_results
    results="test_results/${time}_results.log"
    log="test_results/${time}_test.log"
    coverage="test_results/${time}_coverage.log"
    failed="test_results/${time}_failed_tests_outputs"
    passed="test_results/${time}_passed_tests_outputs"
elif [ $# -eq 5 ]
then
    results=$1
    log=$2
    coverage=$3
    failed=$4
    passed=$5
else
    echo "Wrong number of arguments."
    echo "Usage: $0 <test_results> <test_log> <coverage_log> <failed_test_outputs> <passed_test_outputs>"
    exit 1
fi

which coverage > /dev/null
if [ $? -ne 0 ]; then
    echo "coverage is not installed."
    echo "Install coverage for instance with 'conda install coverage'."
    exit 1
fi

echo "Running python tests."

coverage run --source potential_fitting run_tests.py > $log 2> $results

if [ -d "failed_tests_outputs" ]; then
    echo "Some tests have failed. Outputs of these tests will be found in $failed"
    mv failed_tests_outputs $failed
fi

if [ -d "passed_tests_outputs" ]; then
    echo "The outputs of the tests that passed will be in $passed"
    mv passed_tests_outputs $passed
fi

coverage report -m > $coverage

echo "Done!" 
echo "Test results in ${results}."
echo "Test log in ${log}."
echo "Test covarege report in ${coverage}."
