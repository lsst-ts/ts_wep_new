class Middleware(object):
	
	def __init__(self):
		
		self.filterType = None
		self.pointing = None
		self.cameraRotation = None

	def subscribe(self, item):
		pass

	def publish(self, item, value):
		pass

if __name__ == "__main__":

	pass