#!/usr/bin/env python

import unittest
import StringIO

from cellmlConverter import *


class TestIndenter(unittest.TestCase):
    def setUp(self):
        self.fff = StringIO.StringIO()
        self.out = Indenter(self.fff, "-")
    def tearDown(self):
        self.out = None
        self.fff = None
    def getBuffer(self):
        return self.fff.getvalue()

    def test_simpleOutput(self):
        s = "Testing Testing %s\n"
        self.out(s)
        self.assertEqual(s, self.getBuffer())
    def test_indentCheck(self):
        self.out("0")
        self.out.inc()
        self.out("1")
        self.out("1")
        self.out.inc()
        self.out.dec()
        self.out.inc()
        self.out("2")
        self.out.dec()
        self.out.dec()
        self.out("0")

        self.assertEqual("0\n-1\n-1\n--2\n0\n", self.getBuffer())
    def test_kwargs(self):
        self.out("%(a)s %(b)s %(a)s", b="bbb", a="aaa");
        self.assertEqual("aaa bbb aaa\n", self.getBuffer())
    def test_hasharg(self):
        self.out("%(a)s %(b)s %(a)s", {"a" : "aaa", "b" : "bbb"})
        self.assertEqual("aaa bbb aaa\n", self.getBuffer())
    def test_singleArg(self):
        self.out("xx%sxx", "(A)")
        self.assertEqual("xx(A)xx\n", self.getBuffer())
    def test_multipleArg(self):
        self.out("%s %d %s", "hello", 2, "world")
        self.assertEqual("hello 2 world\n", self.getBuffer())

class TestStandaloneFuncs(unittest.TestCase):
    def test_order(self):
        aaa = ("bbb","ccc","aaa")
        bbb = set(aaa)
        bbb.add("ddd")
        self.assertEqual(tuple(order(bbb)),("aaa","bbb","ccc","ddd"))
        self.assertEqual(tuple(order(aaa)),("aaa","bbb","ccc"))

from xml.etree import cElementTree as ET

def parseEquationXML(sss):
    registerNamespaces()
    eee = stripNamespaces(ET.fromstring(sss))
    return parseEquation(eee)

def printEquation(eqn):
    fff = StringIO.StringIO()
    eqn.toCode(Indenter(fff))
    return fff.getvalue()

class TestEquations(unittest.TestCase):
    def parseEquation(self, sss):
        registerNamespaces()
        eee = stripNamespaces(ET.fromstring(sss))
        return eee
    def printEquation(self, eqn):
        fff = StringIO.StringIO()
        eqn.toCode(Indenter(fff))
        return fff.getvalue()

    def test_istim(self):
        sss = """<apply xmlns="http://www.w3.org/1998/Math/MathML" xmlns:cellml="http://www.cellml.org/cellml/1.0#">
            <eq/>
            <ci>i_Stim</ci>
            <piecewise>
               <piece>
                  <apply>
                     <minus/>
                     <ci>stim_amplitude</ci>
                  </apply>
                  <apply>
                     <and/>
                     <apply>
                        <geq/>
                        <apply>
                           <minus/>
                           <ci>time</ci>
                           <apply>
                              <times/>
                              <apply>
                                 <floor/>
                                 <apply>
                                    <divide/>
                                    <ci>time</ci>
                                    <ci>stim_period</ci>
                                 </apply>
                              </apply>
                              <ci>stim_period</ci>
                           </apply>
                        </apply>
                        <ci>stim_start</ci>
                     </apply>
                     <apply>
                        <leq/>
                        <apply>
                           <minus/>
                           <ci>time</ci>
                           <apply>
                              <times/>
                              <apply>
                                 <floor/>
                                 <apply>
                                    <divide/>
                                    <ci>time</ci>
                                    <ci>stim_period</ci>
                                 </apply>
                              </apply>
                              <ci>stim_period</ci>
                           </apply>
                        </apply>
                        <apply>
                           <plus/>
                           <ci>stim_start</ci>
                           <ci>stim_duration</ci>
                        </apply>
                     </apply>
                  </apply>
               </piece>
               <otherwise>
                  <cn cellml:units="picoA_per_picoF">0</cn>
               </otherwise>
            </piecewise>
         </apply>"""
        self.assertEqual(printEquation(parseEquationXML(sss)), 'i_Stim = (((time-floor(time/stim_period)*stim_period) >= stim_start && (time-floor(time/stim_period)*stim_period) <= (stim_start+stim_duration)) ? -stim_amplitude : 0);\n')

    def test_E_Na(self):
        sss = """         <apply>
            <eq/>
            <ci>E_Na</ci>
            <apply>
               <times/>
               <apply>
                  <divide/>
                  <apply>
                     <times/>
                     <ci>R</ci>
                     <ci>T</ci>
                  </apply>
                  <ci>F</ci>
               </apply>
               <apply>
                  <ln/>
                  <apply>
                     <divide/>
                     <ci>Na_o</ci>
                     <ci>Na_i</ci>
                  </apply>
               </apply>
            </apply>
         </apply>"""        
        self.assertEqual(printEquation(parseEquationXML(sss)), 'E_Na = R*T/F*log(Na_o/Na_i);\n')
    def test_diff_R_prime(self):
        sss = """         <apply xmlns="http://www.w3.org/1998/Math/MathML" xmlns:cellml="http://www.cellml.org/cellml/1.0#">
            <eq/>
            <apply>
               <diff/>
               <bvar>
                  <ci>time</ci>
               </bvar>
               <ci>R_prime</ci>
            </apply>
            <apply>
               <plus/>
               <apply>
                  <times/>
                  <apply>
                     <minus/>
                     <ci>k2</ci>
                  </apply>
                  <ci>Ca_ss</ci>
                  <ci>R_prime</ci>
               </apply>
               <apply>
                  <times/>
                  <ci>k4</ci>
                  <apply>
                     <minus/>
                     <cn cellml:units="dimensionless">1</cn>
                     <ci>R_prime</ci>
                  </apply>
               </apply>
            </apply>
         </apply>"""
        eqn = parseEquationXML(sss);
        eqn.addInitialCondition(Equation("R_prime", "100"));
        self.assertEqual(printEquation(eqn), 'init(R_prime) = 100;\ndiff(R_prime) = (-k2*Ca_ss*R_prime+k4*(1-R_prime));\n')
        #self.assertEqual(self.printEquation(self.parseEquation(sss)), 'E_Na = R*T/F*log(Na_o/Na_i);\n')
        
        
        
