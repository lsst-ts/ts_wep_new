#!/usr/bin/env groovy

pipeline {

    agent {
        // Use the docker to assign the Python version.
        // Use the label to assign the node to run the test.
        // The nodes in T&S teams is 'jenkins-el7-1'.
        // It is recommended by SQUARE team do not add the label.
        docker {
            image 'lsstts/aos:w_2018_47'
            args '-u root'
        }
    }

    triggers {
        pollSCM('H * * * *')
    }

    environment {
        // Use the double quote instead of single quote
        // Add the PYTHONPATH
        PYTHONPATH="${env.WORKSPACE}/python"
        // XML report path
        XML_REPORT="jenkinsReport/report.xml"
        // Module name used in the pytest coverage analysis
        MODULE_NAME="lsst.ts.wep"
    }

    stages {
        stage ('Install Requirements') {
            steps {
                // When using the docker container, we need to change
                // the HOME path to WORKSPACE to have the authority
                // to install the packages.
                withEnv(["HOME=${env.WORKSPACE}"]) {
                    sh """
                        source /opt/rh/devtoolset-6/enable
                        source /opt/lsst/loadLSST.bash
                        conda install scikit-image
                        python builder/setup.py build_ext --build-lib python/lsst/ts/wep/cwfs/lib
                        git clone --branch master https://github.com/lsst/obs_lsst.git
                        cd obs_lsst/
                        git checkout 9c3b73a
                        setup -k -r .
                        scons
                        cd ..
                        git clone --branch master https://github.com/lsst-dm/phosim_utils.git
                        cd phosim_utils/
                        setup -k -r .
                        scons
                    """
                }
            }
        }

        stage('Unit Tests and Coverage Analysis') { 
            steps {
                // Direct the HOME to WORKSPACE for pip to get the
                // installed library.
                // 'PATH' can only be updated in a single shell block.
                // We can not update PATH in 'environment' block.
                // Pytest needs to export the junit report. 
                withEnv(["HOME=${env.WORKSPACE}"]) {
                    sh """
                        source /opt/rh/devtoolset-6/enable
                        source /opt/lsst/loadLSST.bash
                        setup sims_catUtils -t sims_w_2018_47
                        cd obs_lsst/
                        setup -k -r .
                        cd ..
                        cd phosim_utils/
                        setup -k -r .
                        cd ..
                        pytest --cov-report html --cov=${env.MODULE_NAME} --junitxml=${env.WORKSPACE}/${env.XML_REPORT} ${env.WORKSPACE}/tests/cwfs/*.py ${env.WORKSPACE}/tests/bsc/*.py ${env.WORKSPACE}/tests/deblend/*.py ${env.WORKSPACE}/tests/*.py
                    """
                }
            }
        }
    }

    post {        
        always {
            // The path of xml needed by JUnit is relative to
            // the workspace.
            junit 'jenkinsReport/*.xml'

            // Publish the HTML report
            publishHTML (target: [
                allowMissing: false,
                alwaysLinkToLastBuild: false,
                keepAll: true,
                reportDir: 'htmlcov',
                reportFiles: 'index.html',
                reportName: "Coverage Report"
              ])
        }

        cleanup {
            // clean up the workspace
            deleteDir()
        }  
    }
}
