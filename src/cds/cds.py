import pexpect
import skill
import re

class CadenceSession:
	"""Class to handle a non-graphical cadence session

	>>> s = CadenceSession()
	>>> s.send("(list 1 2 3)")
	[1, 2, 3]
	>>> s.callfunc('list', 1,2,3)
	[1, 2, 3]
	>>> s.list(1,2,3)
	[1, 2, 3]
	
	"""
	def __init__(self, cmd="icfb -nograph", timeout=30, verbose=False):
	    self.verbose = verbose
	    self.cds = pexpect.spawn(cmd, timeout=timeout)
	    self.cds.setecho(False)
	    
	    self.prompt = "1> "
	    
	    self.cds.expect(self.prompt)
	    self.startup = self.cds.before
	    if verbose:
		    print self.startup

	def callfunc(self, name, *args, **optargs):
		return self.send("(%s "%name +
				 " ".join(map(skill.toSkill, args) + ["?%s %s"%(k,skill.toSkill(v)) for k,v in optargs.items()])+
				 ")")
	def _makeskillfunc(self, name):
		def func(*args, **optargs):
			return self.callfunc(name, *args, **optargs)
		return func
	
	def __getattr__(self, attr):
		"""Dynamically create a function that executes a skill function"""
		if not attr.startswith("_"):
			return self._makeskillfunc(attr)
		else:
			raise AttributeError()

	def send(self, expr):
		if self.verbose:
			print "Sending: "+expr
		self.cds.sendline(expr)
		self.cds.expect(self.prompt)
		if self.verbose:
			print "Got:", self.cds.before

		response = [s for s in re.split("[\r]\n",self.cds.before) if s != ""][-1]

		try:
		    result = skill.parse(response)
		except:
		    print response
		    raise Exception("Could not parse response")

		if response.startswith("*Error*"):
			raise Exception(response)

		if self.verbose:
			print "Result:", result
		
		return result
	def __del__(self):
		self.cds.sendline("exit")
		self.cds.expect(pexpect.EOF)


if __name__ == "__main__":
    import doctest
    doctest.testmod()

