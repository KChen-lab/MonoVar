import numpy as np
from Single_Cell_Ftrs_Pos import Single_Cell_Ftrs_Pos
from alleles_prior import allele_prior
import math

class Genotype_Prob_matrix:
	def __init__(self, sngl_cell_ftr_list, prior_allele_mat, read_supported_n_cells_others):
		self.n_cells = read_supported_n_cells_others
		self.denom_prob_matrix = np.zeros((2*self.n_cells+1,self.n_cells))
		self.sngl_cell_ftr_list = sngl_cell_ftr_list
		self.prior_allele_mat = prior_allele_mat

        def nCr(self, n, r, factorial_list):
                comb_denom = n-r
                number = factorial_list[n]/(factorial_list[r]*factorial_list[comb_denom])
                return number

	def fill_matrix(self, sngl_cell_ftr_list, original_n_cells, nCr_matrix):
	
		self.denom_prob_matrix[0,0] = float(sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(0))
		
		self.denom_prob_matrix[2,0] = float(sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(2))
		self.denom_prob_matrix[1,0] = 2*float(sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(1))

		for j in range(1,self.n_cells):
			cell_j_prob_0 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(0)
			cell_j_prob_2 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(2)
			cell_j_prob_1 = 2*sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(1)
			for l in range(0, 2*self.n_cells+1):
				if l > 2*(j+1):
					self.denom_prob_matrix[l,j] = 0
				else:
					if (l == 0):
						t1 = 0
						t2 = 0
						t3 = self.denom_prob_matrix[l,j-1]
					elif (l == 1):
						t1 = 0
						t2 = 2*self.denom_prob_matrix[l-1,j-1]
						t3 = self.denom_prob_matrix[l,j-1]
					else:
						t1 = self.denom_prob_matrix[l-2,j-1]
						t2 = 2*self.denom_prob_matrix[l-1,j-1]
                                                t3 = self.denom_prob_matrix[l,j-1]
					# if (sngl_cell_ftr_list[j].depth == 0):
					# 	self.denom_prob_matrix[l,j] = self.denom_prob_matrix[l,j-1] 
					# else:
				#		print sngl_cell_ftr_list[j].prob_reads_given_genotype(prior_allele_mat,2)
				#		print sngl_cell_ftr_list[j].prob_reads_given_genotype(prior_allele_mat,1)
				#		print sngl_cell_ftr_list[j].prob_reads_given_genotype(prior_allele_mat,0)
					self.denom_prob_matrix[l,j] = t1*cell_j_prob_2 + \
								      t2*cell_j_prob_1 + \
						 		      t3*cell_j_prob_0

               	for l in range(0, 2*self.n_cells+1):
                        self.denom_prob_matrix[l,self.n_cells-1] = self.denom_prob_matrix[l,self.n_cells-1]/nCr_matrix[2*self.n_cells, l] 
		
		return self.denom_prob_matrix

	def fill_matrix_stable(self, sngl_cell_ftr_list, original_n_cells, nCr_matrix):
		self.denom_prob_matrix[0,0] = sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(0)
		self.denom_prob_matrix[2,0] = sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(2)
		self.denom_prob_matrix[1,0] = 2*sngl_cell_ftr_list[0].Prob_Reads_Given_Genotype_Genotyping(1)
		sum_l = self.denom_prob_matrix[0,0] + self.denom_prob_matrix[2,0] + self.denom_prob_matrix[1,0]
		self.denom_prob_matrix[0,0] = self.denom_prob_matrix[0,0]/sum_l
		self.denom_prob_matrix[2,0] = self.denom_prob_matrix[2,0]/sum_l
		self.denom_prob_matrix[1,0] = self.denom_prob_matrix[1,0]/sum_l
		for j in range(1,self.n_cells):
			cell_j_prob_0 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(0)
			cell_j_prob_2 = sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(2)
			cell_j_prob_1 = 2*sngl_cell_ftr_list[j].Prob_Reads_Given_Genotype_Genotyping(1)
			sum_l = 0
			for l in range(0, 2*self.n_cells+1):
				
				if (l > 2*(j+1)):
					self.denom_prob_matrix[l,j] = 0
				else:
					if (l == 0):
						t1 = 0
						t2 = 0
						t3 = self.denom_prob_matrix[l,j-1]
					elif (l == 1):
						t1 = 0
						t2 = 2*self.denom_prob_matrix[l-1,j-1]
						t3 = self.denom_prob_matrix[l,j-1]
					else:
						t1 = self.denom_prob_matrix[l-2,j-1]
						t2 = 2*self.denom_prob_matrix[l-1,j-1]
                                                t3 = self.denom_prob_matrix[l,j-1]
					# if (sngl_cell_ftr_list[j].depth == 0):
					# 	self.denom_prob_matrix[l,j] = self.denom_prob_matrix[l,j-1] 
					# else:
					self.denom_prob_matrix[l,j] = t1*cell_j_prob_2 + \
								      t2*cell_j_prob_1 + \
						 		      t3*cell_j_prob_0
					sum_l += self.denom_prob_matrix[l,j]
			for l in range(0, 2*(j+1)):
				self.denom_prob_matrix[l,j] = self.denom_prob_matrix[l,j]/sum_l
					
               	for l in range(0, 2*self.n_cells+1):
                        self.denom_prob_matrix[l,self.n_cells-1] = self.denom_prob_matrix[l,self.n_cells-1]/ nCr_matrix[2*self.n_cells, l] 

        	return self.denom_prob_matrix

#	def printMatrix(self):
#		for l in range (0, 2*self.n_cells+1):
#		#	for j in range (0, n_cells):
#			print self.denom_prob_matrix[l,0], '\t', self.denom_prob_matrix[l,1], '\t', self.denom_prob_matrix[l,2], '\t', \
#			self.denom_prob_matrix[l,3], '\t', self.denom_prob_matrix[l,4], '\t', self.denom_prob_matrix[l,5], '\t', \
#			self.denom_prob_matrix[l,6], '\t', self.denom_prob_matrix[l,7], '\t', self.denom_prob_matrix[l,8], '\t', \
#			self.denom_prob_matrix[l,9], '\t', self.denom_prob_matrix[l,10], '\t', self.denom_prob_matrix[l,11]
