project(metacatalogue)

cmake_minimum_required(VERSION 3.10.0)

add_custom_target(minimap2 COMMAND make 
                   WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/thirdparty/minimap2
)
include(ExternalProject)
ExternalProject_Add(seqwish
        PREFIX seqwish
        SOURCE_DIR ${PROJECT_SOURCE_DIR}/thirdparty/seqwish
        INSTALL_COMMAND "")








