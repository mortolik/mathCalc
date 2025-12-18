cmake_minimum_required(VERSION 3.16)
project(mathCalc LANGUAGES CXX)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt6 COMPONENTS Widgets QUIET)
if(Qt6_FOUND)
    set(QT_PACKAGE Qt6)
else()
    find_package(Qt5 COMPONENTS Widgets REQUIRED)
    set(QT_PACKAGE Qt5)
endif()

add_executable(mathCalc
    src/main.cpp
    src/mainwindow.cpp
    src/mainwindow.h
    src/matrixops.cpp
    src/matrixops.h
)

target_link_libraries(mathCalc PRIVATE ${QT_PACKAGE}::Widgets)
