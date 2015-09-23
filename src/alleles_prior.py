class allele_prior:
	def __init__(self, p):
		self.prior_matrix = {}
		self.prior_matrix[('AA','T')] = p
		self.prior_matrix[('AA','G')] = p
		self.prior_matrix[('AA','C')] = p
		self.prior_matrix[('AA','A')] = 1 - 3*p
		self.prior_matrix[('TT','T')] = 1 - 3*p
	        self.prior_matrix[('TT','G')] = p
	        self.prior_matrix[('TT','C')] = p
	        self.prior_matrix[('TT','A')] = p
                self.prior_matrix[('GG','T')] = p
                self.prior_matrix[('GG','A')] = p
                self.prior_matrix[('GG','C')] = p
                self.prior_matrix[('GG','G')] = 1 - 3*p
                self.prior_matrix[('CC','C')] = 1 - 3*p
                self.prior_matrix[('CC','G')] = p
                self.prior_matrix[('CC','T')] = p
                self.prior_matrix[('CC','A')] = p
		self.prior_matrix[('AC','T')] = p
                self.prior_matrix[('AC','G')] = p
                self.prior_matrix[('AC','C')] = (1 - 2*p)/2
                self.prior_matrix[('AC','A')] = (1 - 2*p)/2
		self.prior_matrix[('AG','T')] = p
                self.prior_matrix[('AG','C')] = p
                self.prior_matrix[('AG','G')] = (1 - 2*p)/2
                self.prior_matrix[('AG','A')] = (1 - 2*p)/2
		self.prior_matrix[('AT','C')] = p
                self.prior_matrix[('AT','G')] = p
                self.prior_matrix[('AT','T')] = (1 - 2*p)/2
                self.prior_matrix[('AT','A')] = (1 - 2*p)/2
                self.prior_matrix[('CG','T')] = p
                self.prior_matrix[('CG','A')] = p
                self.prior_matrix[('CG','C')] = (1 - 2*p)/2
                self.prior_matrix[('CG','G')] = (1 - 2*p)/2
                self.prior_matrix[('CT','A')] = p
                self.prior_matrix[('CT','G')] = p
                self.prior_matrix[('CT','T')] = (1 - 2*p)/2
                self.prior_matrix[('CT','C')] = (1 - 2*p)/2
                self.prior_matrix[('GT','C')] = p
                self.prior_matrix[('GT','A')] = p
                self.prior_matrix[('GT','T')] = (1 - 2*p)/2
                self.prior_matrix[('GT','G')] = (1 - 2*p)/2

	def getValue(self, key):
		return self.prior_matrix[key]

	def PrintPrior(self):
		print self.prior_matrix
