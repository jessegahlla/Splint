import os, re, numpy

#Define exceptions
class TemperatureError(Exception): pass
class InvalidTemperatureInput(TemperatureError): pass
class InvalidExperimentInput(TemperatureError): pass
class InvalidStructureInput(TemperatureError): pass
class InvalidSubmissionInput(TemperatureError): pass
class InvalidExperimentInput(TemperatureError): pass
class InvalidBashInput(TemperatureError): pass
class InvalidWallHours(TemperatureError): pass
class InvalidMDPInput(TemperatureError): pass
class InvalidMDPSettings(TemperatureError): pass

# temperature class
class temperature:
	def __init__(self, name, experiment, structure):
		if not name and len(str(name).split()) == 0:
			raise InvalidTemperatureInput, "Temperature can not be blank"
		if not isinstance(name, float):
			raise InvalidTemperatureInput, "Temperature must be a float"
		if name < 0.0:
			raise InvalidTemperatureInput, "Temperature must be in Kelvins"
		if not experiment:
			raise InvalidExperimentInput, "Experiment can not be blank"
		if not structure:
			raise InvalidStructureInput, "Structure can not be blank"
		if structure.split(".")[-1] != "pdb":
			raise InvalidStructureInput, "Structure must be a .pdb file"
		self.name = name
		self.experiment_name = experiment
		self.structure_file = structure
	def __make_dir__(self):
		os.mkdir(str(self.name))
	def make_mdrun(self, integrator = 'sd', timestep = 0.0005, nsteps = 40000000, nstvout = 0, nstxout = 0, nstlog = 4000, nstenergy = 4000, nstxtcout = 4000, nstlist = 20, rlist = 1.5, rcoulomb = 1.5, rvdw = 1.5, table_extension = 30, tauT = 1.0, refT = 70.0, genTemp = 120, genSeed = 12345):
		#genSeed
		if not genSeed:
			raise InvalidMDPSettings, "genSeed must not be blank"
		if not isinstance(genSeed, int):
			raise InvalidMDPSettings, "genSeed must only be an integer"
		#genTemp
		if not genTemp:
			raise InvalidMDPSettings, "genTemp must not be blank"
		if not (isinstance(genTemp, int) or isinstance(genTemp, float)):
			raise InvalidMDPSettings, "genTemp must only be an integer or float"
		#timestep
		if not timestep:
			raise InvalidMDPSettings, "timestep must not be blank"
		if not (isinstance(timestep, int) or isinstance(timestep, float)):
			raise InvalidMDPSettings, "timestep must only be an integer or float"
		#nstxtcout
		if not nstxtcout:
			raise InvalidMDPSettings, "nstcout must not be blank"
		if not isinstance(nstxtcout, int):
			raise InvalidMDPSettings, "nstcout must only be an integer"
		#nstvout
		if not isinstance(nstvout, int):
			raise InvalidMDPSettings, "nstvout must only be an integer"
		#nstxout
		if not isinstance(nstxout, int):
			raise InvalidMDPSettings, "nstxout must only be an integer"
		#nstlog
		if not nstlog:
			raise InvalidMDPSettings, "nstlog must not be blank"
		if not isinstance(nstlog, int):
			raise InvalidMDPSettings, "nstlog must only be an integer"
		#nstlist
		if not nstlist:
			raise InvalidMDPSettings, "nstlist must not be blank"
		if not isinstance(nstlist, int):
			raise InvalidMDPSettings, "nstlist must only be an integer"
		#nsteps
		if not nsteps:
			raise InvalidMDPSettings, "nsteps must not be blank"
		if not isinstance(nsteps, int):
			raise InvalidMDPSettings, "nsteps must only be an integer"
		#nstenergy
		if not nstenergy:
			raise InvalidMDPSettings, "nstenergy must not be blank"
		if not isinstance(nstenergy, int):
			raise InvalidMDPSettings, "nstenergy must only be an integer"
		#rcoulomb
		if not rcoulomb:
			raise InvalidMDPSettings, "rcoulomb must not be blank"
		if not (isinstance(rcoulomb, int) or isinstance(rcoulomb, float)):
			raise InvalidMDPSettings, "rcoulomb must only be an integer or a float"
		#refT
		if not refT:
			raise InvalidMDPSettings, "refT must not be blank"
		if not (isinstance(refT, int) or isinstance(refT, float)):
			raise InvalidMDPSettings, "refT must only be an integer or a float"
		#tauT
		if not tauT:
			raise InvalidMDPSettings, "tauT must not be blank"
		if not (isinstance(tauT, int) or isinstance(tauT, float)):
			raise InvalidMDPSettings, "tauT must only be an integer or a float"
		#rlist
		if not rlist:
			raise InvalidMDPSettings, "rlist must not be blank"
		if not (isinstance(rlist, int) or isinstance(rlist, float)):
			raise InvalidMDPSettings, "rlist must only be an integer or a float"
		#table-extension
		if not table_extension:
			raise InvalidMDPSettings, "table_extension must not be blank"
		if not (isinstance(table_extension, int) or isinstance(table_extension, float)):
			raise InvalidMDPSettings, "table_extension must only be an integer or a float"
		#rvdw
		if not rvdw:
			raise InvalidMDPSettings, "rvdw must not be a blank"
		if not (isinstance(rvdw, int) or isinstance(rvdw, float)):
			raise InvalidMDPSettings, "rvdw must only be an integer or a float"
		if self.name:
			refT = self.name
		if not os.path.exists(str(refT)):
			os.mkdir(str(refT))
		mdrun = open(str(refT) + "/" + "mdrun.1.mdp", 'a')
		mdrun.write("title\t\t\t= Template\n")
		mdrun.write(";Preprocessor\n")
		mdrun.write("cpp\t\t\t= /lib/cpp\n")
		mdrun.write(";Run control: A leap-from algorithm for integrating Newtonian equations.\n")
		mdrun.write("integrator\t\t\t= %s\n" % integrator)
		mdrun.write(";time step in femtoseconds\n")
		mdrun.write("dt\t\t\t= %s\n" % timestep)
		mdrun.write(";number of stepsi\n")
		mdrun.write("nsteps\t\t\t= %s\n" % nsteps)
		mdrun.write(";frequency to write coordinates to output trajectory file\n")
		mdrun.write("nstxout\t\t\t= %s\n" % nstxout)
		mdrun.write(";frequency to write velocities to output trajectory file\n")
		mdrun.write("nstvout\t\t\t= %s\n" % nstvout)
		mdrun.write(";frequency to write energies to log file\n")
		mdrun.write("nstlog\t\t\t= %s\n" % nstlog)
		mdrun.write(";frequency to write energies to energy file\n")
		mdrun.write("nstenergy\t\t\t= %s\n" % nstenergy)
		mdrun.write(";frequency to write coordinates to xtc trajectory\n")
		mdrun.write("nstxtcout\t\t\t= %s\n" % nstxtcout)
		mdrun.write(";group(s) to write to xtc trajectory\n")
		mdrun.write("xtc_grps\t\t\t= protein\n")
		mdrun.write(";group(s) to write to energy file\n")
		mdrun.write("energygrps\t\t\t= protein\n")
		mdrun.write(";Frequency to update the neighbor list (and the long-range forces,\n;when using twin-range cut-offs).\n")
		mdrun.write("nstlist\t\t\t= %s\n" % nstlist)
		mdrun.write("coulombtype\t\t\t= Cut-off\n")
		mdrun.write(";Make a grid in the box and only check atoms in the neighoring grid cells\n;when constructing a new neighbor list every nstlist steps.\n")
		mdrun.write("ns_type\t\t\t= simple\n")
		mdrun.write(";cut-off distance for the short-range neighbor list\n")
		mdrun.write("rlist\t\t\t= %s\n" % rlist)
		mdrun.write(";treatment of electrostatic interactions\n")
		mdrun.write("rcoulomb\t\t\t= %s\n" % rcoulomb)
		mdrun.write(";treatment of van der waals interactions\n")
		mdrun.write("rvdw\t\t\t= %s\n" % rvdw)
		mdrun.write(";Periodic boundary conditions in all the directions\n")
		mdrun.write("pbc\t\t\t= no\n")
		mdrun.write("table-extension\t\t\t= %s ; (nm) this must be longer than the possible pair interaction\n" % table_extension)
		mdrun.write(";Temperature coupling\n")
		mdrun.write("tc-grps\t\t\t= protein\n")
		mdrun.write("tau_t\t\t\t= %0.1f\n" % tauT)
		mdrun.write("ref_t\t\t\t= %0.2f\n" % refT)
		mdrun.write(";Pressure coupling\n")
		mdrun.write("Pcoupl\t\t\t= no\n")
		mdrun.write(";Velocity generation\n")
		mdrun.write("gen_vel\t\t\t= yes\n")
		mdrun.write("gen_temp\t\t\t= %s\n" % genTemp)
		mdrun.write("gen_seed\t\t\t= %s\n" % genSeed)
		mdrun.write(";Constrain all bonds\n")
		mdrun.write("constraints\t\t\t= none\n")
		mdrun.write("comm_mode\t\t\t= angular\n")
		mdrun.close()
	def make_submission(self, refT = 70.0, structureFile = "structure.pdb", submissionFile = "SUBMIT.1.cmd", experimentName = "experiment", bashFile = "/bin/bash -f", wallHours = 48, mdpFile = "mdrun.1.mdp"):
		if not submissionFile:
			raise InvalidSubmissionInput, "Submission filename must not be blank"
		if not isinstance(submissionFile, str):
			raise InvalidSubmissionInput, "Submission filename must be a string"
		if not experimentName:
			raise InvalidExperimentInput, "Experiment name must not be blank"
		if not isinstance(experimentName, str):
			raise InvalidExperimentInput, "Experiment name must be a string"
		if not isinstance(bashFile, str):
			raise InvalidBashInput, "Bash filename must be a string"
		if not mdpFile:
			raise InvalidMDPInput, "MDP filename must not be blank"
		if not isinstance(mdpFile, str):
			raise InvalidMDPInput, "MDP filename must be a string"
		if not structureFile:
			raise InvalidStructureInput, "Structure name must not be blank"
		if not isinstance(structureFile, str):
			raise InvalidStructureInput, "Structure filename must be a string"
		if structureFile.split(".")[-1] != 'pdb':
			raise InvalidStructureInput, "Structure file must be a .PDB file"
		if not wallHours:
			raise InvalidWallHours, "Wall hours must not be blank"
		if not isinstance(wallHours, int):
			raise InvalidWallHours, "Wall hours must be an integer"
		if self.name and self.structure_file and self.experiment_name:
			refT = self.name
			structureFile = self.structure_file
			experimentName = self.experiment_name
		if not os.path.exists(str(refT)):
			os.mkdir(str(refT))
		submit = open(str(refT) + '/' + submissionFile, 'a')
		refT = str(refT).replace(".", "p")
		submit.write("#!%s\n" % bashFile)
		submit.write("#$ -cwd\n")
		submit.write("#\n")
		submit.write("#$ -N %s.t%s\n" % (self.experiment_name, refT))
		submit.write("#$ -e sge.err\n")
		submit.write("#$ -o sge.out\n")
		submit.write("# requesting %shrs wall clock time\n" % wallHours)
		submit.write("#$ -l h_rt=%s:00:00\n" % wallHours)
		submit.write("#\n")
		submit.write("\n")
		submit.write("\n")
		submit.write("\n")
		submit.write("echo Running on host `hostname`\n")
		submit.write("echo Time is `date`\n")
		submit.write("echo Directory is `pwd`\n")
		submit.write("g_grompp -f %s -c ../%s -p ../%s.top -o %s_t%s.1.tpr\n" % (mdpFile, structureFile, experimentName, experimentName, refT))
		submit.write("g_mdrun -v -s %s_t%s.1.tpr -o %s_t%s_traj.1.trr -c %s_t%s.1.gro -e %s_t%s.1.edr -x %s_t%s.1.xtc\n" % (experimentName, refT, experimentName, refT, experimentName, refT, experimentName, refT, experimentName, refT))

class experiment:
	def __init__(self, name, gro = "structure.pdb"):
		self.name = name
		self.gro = gro
	
	def set_temps(self, temp_list):
		self.temps = {}
		for T in temp_list:
			self.temps[T] = temperature(T, self.name, self.gro)
			
	def build_environ(self):
		for T in self.temps:
			if os.path.exists(str(T)) == False:
				self.temps[T].__make_dir__()
				self.temps[T].make_mdrun()
				self.temps[T].make_submission()