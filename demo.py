import numpy.random as nr
import HC_estimator as hce

def main():
	sample_size = 100
	x_g = nr.uniform(0,1,[sample_size,1])
	y_g = nr.uniform(0,1,[sample_size,1])
	print '*'*100
	print 'bandwidth = 0.53'
	print 'uncorrelated HC:', hce.HC(x_g,y_g,0.53)
	print 'correlated HC', hce.HC(x_g,x_g,0.53)
	print '*'*100
	print 'bandwidth = 1.06'	
	print 'uncorrelated HC:', hce.HC(x_g,y_g,1.06)
	print 'correlated HC:', hce.HC(x_g,x_g,1.06)

if __name__ == '__main__':
	main()
