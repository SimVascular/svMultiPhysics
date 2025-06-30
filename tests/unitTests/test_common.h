#ifndef TEST_COMMON_H
#define TEST_COMMON_H

#include <stdlib.h>
#include <iostream>
#include <random>
#include <chrono>
#include "CepMod.h"
#include "ComMod.h"
#include "gtest/gtest.h"


// --------------------------------------------------------------
// -------------------- Mock svMultiPhysics object -------------------
// --------------------------------------------------------------


class MockCepMod : public CepMod {
public:
    MockCepMod() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockdmnType : public dmnType {
public:
    MockdmnType() {
        // initialize if needed 
    }
    // MockstModelType mockStM;
    // Mock methods if needed
};
class MockmshType : public mshType {
public:
    MockmshType() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockeqType : public eqType {
public:
    MockeqType() {
        // initialize if needed 
    }
    MockdmnType mockDmn;
    // Mock methods if needed
};
class MockComMod : public ComMod {
public:
    MockComMod() {
        // initialize if needed 
        nsd = 3;
    }
    MockeqType mockEq;
    MockmshType mockMsh;
    // Mock methods if needed
};


// --------------------------------------------------------------
// --------------------Base class for testing ------------------------
// --------------------------------------------------------------
class TestBase {
public:
    MockComMod com_mod;
    MockCepMod cep_mod;
};


#endif