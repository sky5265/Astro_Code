import pandas as pd


def get_element_symbol(atomic_number):
	df = pd.read_csv('Periodic_Table.csv')

	element = df.loc[df['AtomicNumber'] == atomic_number]
	l = list(element['Symbol'])
	if len(l) == 1:
		return l[0]
	else:
		raise ValueError("Atomic Number: "+str(atomic_number)+" is invalid!")
