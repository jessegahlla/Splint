#!/usr/bin/env python
# encoding: utf-8
"""
splint_temperature_test.py

Created by Jesse Singh on 2011-03-23.
Copyright (c) 2011 Jesse Singh. All rights reserved.

temperature should return an object under correct conditions
temperature.__make_dir__ should create a directory under correct conditions
temperature.make_mdrun should create an mdrun file under correct conditions
temperature.make_submission should create a submission file under correct conditions

temperature should fail if name is not a number
temperature should fail if name is not a float
temperature should fail if name is less than 0
temperature should fail if name is blank
temperature should fail if experiment name is blank
temperature should fail if gro is not a .gro or .pdb file
temperature should fail if gro is blank

temperature.make_mdrun should fail if "name" directory does not exist

temperature.make_submission should fail if "submission" directory does not exist
temperature.make_submission should fail if structureFile is not a string
temperature.make_submission should fail if structureFile does not end in pdb
temperature.make_submission should fail if structureFile is blank
temperature.make_submission should fail if submissionFile is not a string
temperature.make_submission should fail if submissionFile does not end in cmd
temperature.make_submission should fail if submissionFile is blank
temperature.make_submission should fail if experimentName is not a string
temperature.make_submission should fail if experimentName is blank
temperature.make_submission should fail if bashFile is not a string
temperature.make_submission should fail if wallHours is not an integer
temperature.make_submission should fail if wallHours is blank
temperature.make_submission should fail if mdpFile is not a string
temperature.make_submission should fail if mdpFile is blank
"""
import splint
import os
import unittest

class TemperatureKnownValues(unittest.TestCase):
	knownValues = ( (0.0, "mixZZ100100", "mixZZ100100.pdb"),
				(1.0, "conversion", "structure.pdb"),
				(5.23, "test_case5860", "test_5860.pdb"),
				(10.7, "mix.test", "mix.pdb"),
				(100.76, "conversionCase", "structure.pdb"),
				(117.996, "ZZ125100", "protein.pdb"),
				(356.3529, "ab100100", "ab100100.pdb"))
	
	knownMDRun = ["title\t\t\t= Template\n",
				  ";Preprocessor\n",
				  "cpp\t\t\t= /lib/cpp\n",
				  ";Run control: A leap-from algorithm for integrating Newtonian equations.\n",
				  "integrator\t\t\t= sd\n",
				  ";time step in femtoseconds\n",
				  "dt\t\t\t= 0.0005\n",
				  ";number of stepsi\n",
				  "nsteps\t\t\t= 40000000\n",
				  ";frequency to write coordinates to output trajectory file\n",
				  "nstxout\t\t\t= 0\n",
				  ";frequency to write velocities to output trajectory file\n",
				  "nstvout\t\t\t= 0\n",
				  ";frequency to write energies to log file\n",
				  "nstlog\t\t\t= 4000\n",
				  ";frequency to write energies to energy file\n",
				  "nstenergy\t\t\t= 4000\n",
				  ";frequency to write coordinates to xtc trajectory\n",
				  "nstxtcout\t\t\t= 4000\n",
				  ";group(s) to write to xtc trajectory\n",
				  "xtc_grps\t\t\t= protein\n",
				  ";group(s) to write to energy file\n",
				  "energygrps\t\t\t= protein\n",
				  ";Frequency to update the neighbor list (and the long-range forces,\n",
				  ";when using twin-range cut-offs).\n",
				  "nstlist\t\t\t= 20\n",
				  "coulombtype\t\t\t= Cut-off\n",
				  ";Make a grid in the box and only check atoms in the neighoring grid cells\n",
				  ";when constructing a new neighbor list every nstlist steps.\n",
				  "ns_type\t\t\t= simple\n",
				  ";cut-off distance for the short-range neighbor list\n",
				  "rlist\t\t\t= 1.5\n",
				  ";treatment of electrostatic interactions\n",
				  "rcoulomb\t\t\t= 1.5\n",
				  ";treatment of van der waals interactions\n",
				  "rvdw\t\t\t= 1.5\n",
				  ";Periodic boundary conditions in all the directions\n",
				  "pbc\t\t\t= no\n",
				  "table-extension\t\t\t= 30 ; (nm) this must be longer than the possible pair interaction\n",
				  ";Temperature coupling\n",
				  "tc-grps\t\t\t= protein\n",
				  "tau_t\t\t\t= 1.0\n",
				  "ref_t\t\t\t= 70.00\n",
				  ";Pressure coupling\n",
				  "Pcoupl\t\t\t= no\n",
				  ";Velocity generation\n",
				  "gen_vel\t\t\t= yes\n",
				  "gen_temp\t\t\t= 120\n",
				  "gen_seed\t\t\t= 12345\n",
				  ";Constrain all bonds\n",
				  "constraints\t\t\t= none\n",
				  "comm_mode\t\t\t= angular\n"]

	knownSubmission = ["#!/bin/bash -f\n",
					   "#$ -cwd\n",
					   "#\n",
					   "#$ -N test.t70p2\n",
					   "#$ -e sge.err\n",
					   "#$ -o sge.out\n",
					   "# requesting 48hrs wall clock time\n",
					   "#$ -l h_rt=48:00:00\n",
					   "#\n",
					   "\n",
					   "\n",
					   "\n",
					   "echo Running on host `hostname`\n",
					   "echo Time is `date`\n",
					   "echo Directory is `pwd`\n",
					   "g_grompp -f mdrun.1.mdp -c ../structure.pdb -p ../test.top -o test_t70p2.1.tpr\n",
					   "g_mdrun -v -s test_t70p2.1.tpr -o test_t70p2_traj.1.trr -c test_t70p2.1.gro -e test_t70p2.1.edr -x test_t70p2.1.xtc\n"]

	def testTemperatureKnownValues(self):
		"""temperature objects should be created with known input"""
		for temp, experiment, structure in self.knownValues:
			temperatureTestObject = splint.temperature(temp, experiment, structure)
			self.assertEqual(temperatureTestObject.name, temp)
			self.assertEqual(temperatureTestObject.experiment_name, experiment)
			self.assertEqual(temperatureTestObject.structure_file, structure)
	
	def testTemperatureMakeDir(self):
		"""temperature directories should be created with known input"""
		for temp, experiment, structure in self.knownValues:
			temperatureTestObject = splint.temperature(temp, experiment, structure)
			if not os.path.isdir(str(temp)):
				temperatureTestObject.__make_dir__()
				tempDir = str(temp)
				self.assertRaises(OSError, os.mkdir, tempDir)
				os.rmdir(str(temp))
			
	def testTemperatureMakeMDRUN(self):
		"""mdrun files should be created with known input"""
		testTemp = 70.0
		temperatureTestObject = splint.temperature(testTemp, "test_experiment", "structure.pdb")
		if not os.path.isdir(str(testTemp)):
			os.mkdir(str(testTemp))
			temperatureTestObject.make_mdrun()
			f = open(str(testTemp) + '/' + 'mdrun.1.mdp', 'r')
			i = 0
			for line in f:
				self.assertEqual(line, self.knownMDRun[i])
				i += 1
			f.close()
			os.remove(str(testTemp) + '/' + 'mdrun.1.mdp')
			os.rmdir(str(testTemp))

	def testTemperatureMakeSubmission(self):
		"""Submission files should be created with known input"""
		testTemp = 70.2
		temperatureTestObject = splint.temperature(testTemp, "test", "structure.pdb")
		if not os.path.isdir(str(testTemp)):
			os.mkdir(str(testTemp))
			temperatureTestObject.make_submission()
			f = open(str(testTemp) + '/' + 'SUBMIT.1.cmd', 'r')
			i = 0
			for line in f:
				self.assertEqual(line, self.knownSubmission[i])
				i += 1
			f.close()
			os.remove(str(testTemp) + '/' + 'SUBMIT.1.cmd')
			os.rmdir(str(testTemp))
		
class TemperatureBadInput(unittest.TestCase):
		
	def testNameNotNumber(self):
		"""temperature should fail if name is not a number"""
		self.assertRaises(splint.InvalidTemperatureInput, splint.temperature, "seventy", "mixZZ100100", "mixZZ100100.pdb")

	def testNameNotFloat(self):
		"""temperature should fail if name is not a float"""
   		self.assertRaises(splint.InvalidTemperatureInput, splint.temperature, 70, "mixZZ100100", "mixZZ100100.pdb")

	def testNameZero(self):
		"""temperature should fail if name is less than 0.0"""
		self.assertRaises(splint.InvalidTemperatureInput, splint.temperature, -1.0, "mixZZ100100", "mixZZ100100.pdb")

	def testNameBlank(self):
		"""temperature should fail if no name is given"""
		self.assertRaises(splint.InvalidTemperatureInput, splint.temperature, "", "mixZZ100100", "mixZZ100100.pdb")

	def testExperimentBlank(self):
		"""temperature should fail if experiment name is blank"""
		self.assertRaises(splint.InvalidExperimentInput, splint.temperature, 70.0, "", "mixZZ100100.pdb")

	def testStructureNotPDB(self):
		"""temperature should fail if structurefile is not a .pdb file"""
		self.assertRaises(splint.InvalidStructureInput, splint.temperature, 70.0, "mixZZ100100", "mixZZ100100.gro")

	def testStructureBlank(self):
		"""temperature should fail if structure name is blank"""
		self.assertRaises(splint.InvalidStructureInput, splint.temperature, 70.0, "mixZZ100100", "")		

class TemperatureMakeSubmissionBadInput(unittest.TestCase):
	
	def setUp(self):
		testTemperature = 0.1
		testExperiment = "testObject"
		testStructure = "testStructure.pdb"
		self.testTempObject = splint.temperature(testTemperature, testExperiment, testStructure)

	def tearDown(self):
		if os.path.isdir(str(0.1)):
			if os.path.exists(str(0.1) + '/' + 'SUBMIT.1.cmd'): os.remove(str(0.1) + '/' + 'SUBMIT.1.cmd')
			os.rmdir(str(0.1))
		
	def testStructureFileNotString(self):
		"""temperature.make_submission should fail if structureFile is not a string"""
		self.assertRaises(splint.InvalidStructureInput, self.testTempObject.make_submission, structureFile = 0.4)

	def testStructureFileNotPDB(self):
		"""temperature.make_submission should fail if structureFile does not end in pdb"""
		self.assertRaises(splint.InvalidStructureInput, self.testTempObject.make_submission, structureFile = "structure.gro")
	
	def testStructureFileBlank(self):
		"""temperature.make_submission should fail if structureFile is blank"""
		self.assertRaises(splint.InvalidStructureInput, self.testTempObject.make_submission, structureFile = "")
		
	def testSubmissionFileNotString(self):
		"""temperature.make_submission should fail if submissionFile is not a string"""
		self.assertRaises(splint.InvalidSubmissionInput, self.testTempObject.make_submission, submissionFile = 0.4)

	def testSubmissionFileBlank(self):
		"""temperature.make_submission should fail if submissionFile is blank"""
		self.assertRaises(splint.InvalidSubmissionInput, self.testTempObject.make_submission, submissionFile = "")

	def testExperimentNotString(self):
		"""temperature.make_submission should fail if experimentName is not a string"""
		self.assertRaises(splint.InvalidExperimentInput, self.testTempObject.make_submission, experimentName = 0.7)
	
	def testExperimentBlank(self):
		"""temperature.make_submission should fail if experimentName is blank"""
		self.assertRaises(splint.InvalidExperimentInput, self.testTempObject.make_submission, experimentName = "")
	
	def testBashFileNotString(self):
		"""temperature.make_submission should fail if bashFile is not a string"""
		self.assertRaises(splint.InvalidBashInput, self.testTempObject.make_submission, bashFile = 0.7)

	def testWallHoursNotInteger(self):
		"""temperature.make_submission should fail if wallHours is not an integer"""
		self.assertRaises(splint.InvalidWallHours, self.testTempObject.make_submission, wallHours = "twelve")
	
	def testWallHoursNotBlank(self):
		"""temperature.make_submission should fail if wallHours is blank"""
		self.assertRaises(splint.InvalidWallHours, self.testTempObject.make_submission, wallHours = "")

	def testMDPFileNotString(self):
		"""temperature.make_submission should fail if mdpFile is not a string"""
		self.assertRaises(splint.InvalidMDPInput, self.testTempObject.make_submission, mdpFile = 1234)

	def testMDPFileBlank(self):		
		"""temperature.make_submission should fail if mdpFile is blank"""
		self.assertRaises(splint.InvalidMDPInput, self.testTempObject.make_submission, mdpFile = "")

class TemperatureMakeMDRunBadInput(unittest.TestCase):

	def setUp(self):
		testTemperature = 0.1
		testExperiment = "testObject"
		testStructure = "testStructure.pdb"
		self.testTempObject = splint.temperature(testTemperature, testExperiment, testStructure)

	def tearDown(self):
		if os.path.isdir(str(0.1)):
			if os.path.exists(str(0.1) + '/' + 'mdrun.1.mdp'): os.remove(str(0.1) + '/' + 'mdrun.1.mdp')
			os.rmdir(str(0.1))

	def testTimestepNotFloat(self):
		"""temperature.make_mdrun should fail if timestep is not an integer or a float"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, timestep = "0.003")

	def testTimestepBlank(self):
		"""temperature.make_mdrun should fail if timestep is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, timestep = "")
		
	def testNstepsNotInt(self):
		"""temperature.make_mdrun should fail if nsteps is not a int"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nsteps = 400.23)

	def testNstepsBlank(self):
		"""temperature.make_mdrun should fail if nsteps is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nsteps = "")

	def testNstxoutNotInt(self):
		"""temperature.make_mdrun should fail if nstxout is not a int"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstxout = 400.23)

	def testNstxoutBlank(self):
		"""temperature.make_mdrun should fail if nstxout is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstxout = "")

	def testNstvoutNotInt(self):
		"""temperature.make_mdrun should fail if nstvout is not a int"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstvout = 400.23)

	def testNstvoutBlank(self):
		"""temperature.make_mdrun should fail if nstvout is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstvout = "")

	def testNstlogNotInt(self):
		"""temperature.make_mdrun should fail if nstlog is not a int"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstlog = 400.23)

	def testNstlogBlank(self):
		"""temperature.make_mdrun should fail if nstlog is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstlog = "")

	def testNstenergyNotInt(self):
		"""temperature.make_mdrun should fail if nstenergy is not a int"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstenergy = 400.23)

	def testNstenergyBlank(self):
		"""temperature.make_mdrun should fail if nstenergy is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstenergy = "")

	def testNstcoutNotInt(self):
		"""temperature.make_mdrun should fail if nstcout is not a int"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstxtcout = 400.23)

	def testNstcoutBlank(self):
		"""temperature.make_mdrun should fail if nstcout is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstxtcout = "")
		
	def testNstlistNotInt(self):
		"""temperature.make_mdrun should fail if nstlist is not a int"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstlist = 400.23)

	def testNstlistBlank(self):
		"""temperature.make_mdrun should fail if nstlist is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, nstlist = "")

	def testRlistString(self):
		"""temperature.make_mdrun should fail if rlist is not an int or float"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, rlist = "twenty")

	def testRlistBlank(self):
		"""temperature.make_mdrun should fail if rlist is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, rlist = "")

	def testRcoulombString(self):
		"""temperature.make_mdrun should fail if rcoulomb is not an int or float"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, rcoulomb = "twenty")

	def testRcoulombBlank(self):
		"""temperature.make_mdrun should fail if rcoulomb is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, rcoulomb = "")

	def testRvdwString(self):
		"""temperature.make_mdrun should fail if rvdw is not an int or float"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, rvdw = "twenty")

	def testRvdwBlank(self):
		"""temperature.make_mdrun should fail if rvdw is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, rvdw = "")

	def testTableExtensionString(self):
		"""temperature.make_mdrun should fail if table_extension is not an int or float"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, table_extension = "twenty")

	def testTableExtensionBlank(self):
		"""temperature.make_mdrun should fail if table_extension is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, table_extension = "")

	def testTauTString(self):
		"""temperature.make_mdrun should fail if tauT is not an int or float"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, tauT = "twenty")

	def testTauTBlank(self):
		"""temperature.make_mdrun should fail if tauT is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, tauT = "")

	def testRefTString(self):
		"""temperature.make_mdrun should fail if refT is not an int or float"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, refT = "twenty")

	def testGenTempString(self):
		"""temperature.make_mdrun should fail if genTemp is not an int or float"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, genTemp = "twenty")

	def testGenTempBlank(self):
		"""temperature.make_mdrun should fail if genTemp is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, genTemp = "")

	def testGenSeedString(self):
		"""temperature.make_mdrun should fail if genSeed is a string"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, genSeed = "twenty")

	def testGenSeedFloat(self):
		"""temperature.make_mdrun should fail if genSeed is a float"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, genSeed = 213.12)

	def testGenSeedBlank(self):
		"""temperature.make_mdrun should fail if genSeed is blank"""
		self.assertRaises(splint.InvalidMDPSettings, self.testTempObject.make_mdrun, genSeed = "")

if __name__ == '__main__':
	unittest.main()