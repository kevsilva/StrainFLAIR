project(metacatalogue)

cmake_minimum_required(VERSION 3.10.0)


include(ExternalProject)
ExternalProject_Add(minimap2
        PREFIX minimap2
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/minimap2
        INSTALL_COMMAND "")





