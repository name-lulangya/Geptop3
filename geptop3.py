#!/usr/bin/python
# -*- coding: utf-8 -*-
import os
import sys
import math  
from Bio import SeqIO
from Bio.Blast import NCBIXML 
import pickle	
from multiprocessing import Pool
import re
import random
import numpy as np
import pandas as pd
from fractions import Fraction
from imblearn.over_sampling import SMOTE
from sklearn.svm import SVC
from sklearn.model_selection import cross_val_score, GridSearchCV
from sklearn.linear_model import LinearRegression
import joblib
from Bio.Seq import translate


def getargu():
	"""Get argument from terminal"""
	LenArgv = len(sys.argv) 
	Cutoff = 0.3 # default value
	for i in range(LenArgv):
		if sys.argv[i] == '-i':
			Target=sys.argv[i+1]  
		if sys.argv[i] == '-s':
			Cutoff=float(sys.argv[i+1])   
	return Target, Cutoff

def translation(Target, userFileName):
	"""Translate nucleic acid sequence files to protein sequence files"""
	outFile = "{}/translatedFile/{}".format(rootdir,userFileName)
	nucl_dict = SeqIO.to_dict(SeqIO.parse(Target,"fasta"), key_function = lambda rec: rec.description) # keep whitespace in FASTA header
	prot_dict = {}
	num = 0
	for key in nucl_dict:
		nucl_seq = nucl_dict[key].seq
		if len(nucl_seq)%3==0:
			prot_seq = translate(nucl_seq,stop_symbol="")
			prot_dict[key] = prot_seq
		else:
			num+=1
	print("The number of genes whose number of bases is not a multiple of three is {}".format(num))
	with open("{}".format(outFile),'w') as result:
		for k,v in prot_dict.items():
			s = re.sub(r"(.{70})", "\\1\r\n", str(v))
			result.write(">"+"{}".format(k)+"\r\n"+"{}".format(s)+"\r\n")
	return outFile

def information_exstract(Target_Record):
	"""Extraction of Essential Information""" 
	Geptop_prediction = {}  # Store the scores of homologous mapping
	GeneID = []  # Store gene id

	for record in Target_Record:
		if "lcl|" in record.id:
			Geptop_prediction[record.id.split("|")[1]] = 0
			GeneID.append(record.id.split("|")[1])
		else:
			Geptop_prediction[record.id] = 0
			GeneID.append(record.id)
	
	deg=open('{}/DEG3'.format(rootdir),'r')
	EssentialGeneid=[]  # Gene id od reference species
	for row in deg:
		if "|"in row:
			EssentialGeneid.append(row.split("|")[1])  
		else:
			EssentialGeneid.append(row.strip())   
	deg.close()

	return Geptop_prediction, GeneID, EssentialGeneid

def StrToNum(string):
	aa_dict={'A':0, 'C':1, 'D':2, 'E':3, 'F':4, 'G':5, \
			 'H':6, 'I':7, 'K':8, 'L':9, 'M':10, 'N':11, \
			 'P':12, 'Q':13, 'R':14, 'S':15, 'T':16, 'V':17, \
			 'W':18, 'Y':19, 'B':2, 'U':1, 'X':5, 'Z':3,'J':7}	
	number=0
	for i in range(len(string)):
		if string[i] not in aa_dict:
			return -1
		bit=aa_dict[string[i]]
		power=len(string)-i-1
		number+=bit*(20**power)
	return number

def CompositionVector(genomeSeq):
	"""Get the composition vector of species"""
	kstring=6
	k={}
	k0={}
	k1={}
	k2={}
	try:
		for seqrecord in SeqIO.parse(genomeSeq,"fasta"):
			seq=str(seqrecord.seq)
			len0=len(seq)-kstring+3
			for s in range(len0):
				start=s
				end=kstring+s-2
				num=StrToNum(seq[start:end])
				if num not in k2:
					k2[num]=1
				else: k2[num]+=1
				if s<len0-2:
					num=StrToNum(seq[start:end+2])
					if num not in k0:
						k0[num]=1
					else: k0[num]+=1
				if s<len0-1:
					num=StrToNum(seq[start:end+1])
					if num not in k1:
						k1[num]=1
					else: k1[num]+=1
			if -1 in k0:del k0[-1]
			if -1 in k1:del k1[-1]
			if -1 in k2:del k2[-1]
			string0=sum(k0.values())
			string1=sum(k1.values())
			string2=sum(k2.values())
		for n1 in k1:
			for aa in range(20):
				n0=n1*20+aa
				n2=n1%(20**(kstring-2))*20+aa
				if n2 in k1:
					if n0 in k0:
						n3=n1%(20**(kstring-2))
						p0=1.0*k1[n1]*k1[n2]*string0*string2/k2[n3]/string1/string1
						k[n0]=(k0[n0]-p0)/p0
					else:k[n0]=-1
		return k
	except IOError as err:
			print('Compute CV failed to: ',err)

def Distance(CV1,CV2):
	"""Calculate the similarity between two species based on the combined composition"""
	O=0
	P=0
	Q=0
	for value in CV1:
		if value in CV2:
			O+=CV1[value]*CV2[value]
			P+=CV1[value]*CV1[value]
			Q+=CV2[value]*CV2[value]
		else:P+=CV1[value]*CV1[value]
	for value in CV2:
		if value not in CV1:
			Q+=CV2[value]*CV2[value]
	return O/math.sqrt(P*Q)

def BLAST(Target,OGfaa):
	"""Obtain orthologous genes through BLASTP"""
	Tfilename=os.path.split(Target)[-1]
	Ofilename=os.path.split(OGfaa)[-1]
	
	outxml=os.path.join("{}/xmlresult".format(rootdir),Tfilename+'_'+Ofilename+'.xml')
	try:
		BLAST_Results=[]

		blastp='blastp -query '+Target+' -db '+os.path.join("{}/blastdatabase".format(rootdir),Ofilename)+' -out '+outxml+' -outfmt 5 -num_threads 8'
		print(blastp)
		os.system(blastp)
		
		Blast_Records=NCBIXML.parse(open(outxml))
		for blast_record in Blast_Records:
			if "lcl|" in blast_record.query:
				Query_ID=blast_record.query.split("|")[1].split(" ")[0]
			else:
				Query_ID=blast_record.query.split(" ")[0]
			for alignment in blast_record.alignments:
				Hit_ID=alignment.hit_id
				
				for hsp in alignment.hsps:
					E_value=hsp.expect
					if E_value<10:
						BLAST_Results.append((str(Query_ID),str(Hit_ID)))
						break
				break
	except Exception as err:
		print('BLAST failed to: ',err)
	finally: 
		os.remove(outxml)
		return (BLAST_Results)
	
def homology_mapping(TargetCV,Target,OGfaa):
	"""Workflow of homology mapping"""
	blast_result=[]
	with open(os.path.join('{}/CVFile'.format(rootdir),os.path.split(OGfaa)[-1]+".pkl"), "rb") as f:
		DataSetCV = pickle.load(f)	
	SpeciesDistance=Distance(TargetCV,DataSetCV)

	BLAST_Results1=BLAST(Target,OGfaa)
	BLAST_Results2=BLAST(OGfaa,Target)

	for (query,hit) in BLAST_Results1:
		if (hit,query) in BLAST_Results2:
			blast_result.append((query,hit))
		elif (hit,f"ref|{query}|") in BLAST_Results2:
			blast_result.append((query,hit))
		elif (hit.split("|")[1],query) in BLAST_Results2:
			blast_result.append((query,hit))
		elif (hit.split("|")[1],f"ref|{query}|") in BLAST_Results2:
			blast_result.append((query,hit))

	return (OGfaa,blast_result,SpeciesDistance)

def homo_result_extract(result, Geptop_prediction, EssentialGeneid):
	"""Obtain the results of homology mapping"""
	dataset = ""
	blast_result = []
	distance_result = 0
	spieces_distance_result = []   

	for work in result:
		(dataset,blast_result,distance_result) = work.get()
		spieces_distance_result.append([dataset,distance_result])
		for (query,hit) in blast_result:
			if hit.split('|')[1] in EssentialGeneid:
				a = 1.0*distance_result
				Geptop_prediction[query]=Geptop_prediction[query]+a
		print(dataset + "_has been done!")

	return spieces_distance_result

def result1(Geptop_prediction, GeneID):
	"""Return the normalized results of homology mapping"""
	Geptop_max_predict = max(Geptop_prediction.values())
	Geptop_min_predict = min(Geptop_prediction.values())

	print(Geptop_max_predict,Geptop_min_predict)

	Geptop_score = {}  
	for pid in GeneID:
		if Geptop_max_predict == Geptop_min_predict:
			geptop_score = Geptop_prediction[pid]
		else:
			geptop_score = (Geptop_prediction[pid] - Geptop_min_predict) / (Geptop_max_predict-Geptop_min_predict)   
		Geptop_score[pid] = geptop_score  

	return Geptop_score

def remove_file(userFileName):
	"""Delete the created library files"""
	suffix=['.pdb','.phd','.phi','.phr','.pin','.pnd','.pni','.pog','.pos','.pot','.psq','.ptf','.pto']
	for hz in suffix:
		rmfile=os.path.join('{}/blastdatabase'.format(rootdir),userFileName+hz)
		if os.path.exists(rmfile):
			os.remove(rmfile)
	print("remove_file done!")

def seed_select(Geptop_result):
	"""Select the seeds of the sequence composition model"""
	geptop_sorted = sorted(Geptop_result.items(), key = lambda e:e[1],reverse=True)
	geptop_id = [g[0] for g in geptop_sorted]
	geptop_e_num = 120
	geptop_ne_num = int(-0.2*len(geptop_id))
	id_essential = geptop_id[:geptop_e_num]
	id_nonessential = geptop_id[geptop_ne_num:]
	print("seed_select done!")

	return id_essential, id_nonessential

def id_seq_get(userFileName):
	"""Obtain the id-seq dictionary of nucleic acid FASTA file"""
	id_seq = {}
	nuc_handle = open("{}/uploadFile/{}".format(rootdir, userFileName),"r")
	for c in SeqIO.parse(nuc_handle,"fasta"):
		if "lcl|" in c.id:
			id_seq[c.id.split("|")[1]] = c.seq.upper()
		else:
			id_seq[c.id] = c.seq.upper()
	return id_seq

def multi_process_feature_get(id_seq, n):
	"""Extract sequence features using n processes"""
	length = len(id_seq)
	id_seq_list = list(id_seq.items())
	b = [dict(id_seq_list[i:i+n]) for i in range(0, length, n)]
	p=Pool(n)  # Number of parallel processes
	result = []
	for each_id_seq in b:
		result.append(p.apply_async(feature_get,args=(each_id_seq,)))
	p.close()
	p.join()
	feature_list = []
	for work in result:
		each_feature, each_delete = work.get()
		feature_list.append(each_feature)
	merged_df = pd.concat(feature_list, ignore_index=False)
	return merged_df
def feature_get(id_seq):
	"""Extraction of sequence frequency feature"""
	K, L = 4,3
	N = 3*(4+16+64+256)
	delete = []
	col = np.arange(N*L-4*L*L+N)
	df = pd.DataFrame(columns=col)
	df.index.name = "id"
	for k,v in id_seq.items():
		data = feature_extract(k, v, K, L, N, delete)
		df.loc[k] = data
	return df, delete

def feature_extract(id,seq, K, L, N, delete):
	feature = [0]*(N*L-4*L*L+N+3)
	seq_list = list(seq)
	base_num = {'A':1, 'T':2, 'C':3, 'G':4, 'R':random.choice([1,4]), 'Y':random.choice([2,3]), 'M':random.choice([1,3]),\
				'K':random.choice([2,4]), 'S':random.choice([3,4]), 'W':random.choice([1,2]), 'H':random.choice([1,2,3]),\
					'B':random.choice([2,3,4]), 'V':random.choice([1,3,4]), 'D':random.choice([1,2,4]), 'N':random.choice([1,2,3,4])}
	G = len(seq_list)
	try:
		if G%3!=0:
			delete.append(id)
			raise ValueError("Not triplets !")
		else:
			for k in np.arange(1, K+1):
				for l in np.arange(0, L+1):
					if k == 1 and l > 0:
						break
					else:
						for s in np.arange(G):
							if s+k+l <= (G):
								if (k-1+l)%3 == 1 and s%3== 0 or 1: 
									d = int((G-k+1-l)/3+1)
								elif (k-1+l)%3 == 2 and s%3== 0: 
									d = int((G-k+1-l)/3+1)
								else:
									d = int((G-k+1-l)/3)
								#print(d)
								if d !=0:
									feature[caculate(k,l,s,seq_list,base_num,N)] = feature[caculate(k,l,s,seq_list,base_num,N)]+Fraction(1,d)
								else:
									break
							else:
								break
	except ValueError as v:
		print("{}{}".format(id, v))
	del feature[:3]
	feature = [float(x) for x in feature]
	return feature

def caculate(k,l,s,seq_list,base_num,N):
	num = [0]*(k)
	for n in np.arange(0,k-1):
		num[n] = base_num[seq_list[s+n]]
	num[k-1] = base_num[seq_list[s+l+k-1]]
	a = 0
	for m in np.arange(len(num)):
		a += num[m]*4**m
	a_finally = 3*a+l*N-4*l*l+s%3
	return a_finally

def dataframe_drop(df, delete):
	"""Remove the genes whose number of bases is not a multiple of three"""
	df1 = df.drop(index = delete)
	return df1

def seed_feature_extract(df, seed_essential, seed_nonessential):
	"""Add the essentiality label to the seeds"""
	es = []
	seed_id = seed_essential + seed_nonessential
	print("len(seed_id):{}".format(len(seed_id)))
	seed_df = df.loc[seed_id]
	seed_df_id = list(seed_df.index)
	for GI in seed_df_id:
		if GI in seed_essential:
			es.append(1)
		else:
			es.append(0)
	seed_df.insert(len(seed_df.columns),"essensiality",es)
	return seed_df

def train_model(X_ban, Y_ban):
    """Train the Support Vector Machine model using grid search"""
    prammar_dict = {'C': [0.1,0.5,1,5,10,50,100,500,1000], 'gamma': [0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5,1,5,10],
                    'kernel': ['rbf'], 'class_weight':['balanced']}
    svc = SVC(probability=True, random_state= 306)
    grid = GridSearchCV(estimator = svc ,param_grid = prammar_dict, cv = 5, refit = True,scoring = 'roc_auc', n_jobs = 10,)
    result = grid.fit( X_ban, Y_ban)

    print('result.best_score:',result.best_score_)  
    print('result.best_params:',result.best_params_)  
    print('result.best_estimator:',result.best_estimator_)  
    print('model trained')
    return result.best_estimator_

def prediction(model, feature):
    """Make predictions using the sequence composition model"""
    scfm_id_score = {}
    data = feature
    input_data = data.iloc[:,1:4044]
    gene_id = list(data.index)
    score_predict = model.predict_proba(input_data)
    score_result = score_predict[:,1]
    for f in range(len(gene_id)):
        scfm_id_score[gene_id[f]] = score_result[f]
    return scfm_id_score

def coefficient_caculate(userFileName, seed_5folds_auc, spieces_distance_result):
	"""Calculate the coefficients of the two models"""
	geptop_model = LinearRegression()
	geptop_params = np.load("{}/geptop_auc_predict_LR.npz".format(rootdir))
	geptop_model.coef_ = geptop_params['coef']
	geptop_model.intercept_ = geptop_params['intercept']
	
	scfm_model = LinearRegression()
	scfm_params = np.load("{}/scfm_auc_predict_LR.npz".format(rootdir))
	scfm_model.coef_ = scfm_params['coef']
	scfm_model.intercept_ = scfm_params['intercept']

	geptop_aucTure_file = open("{}/reference_auc.txt".format(rootdir),'r')
	geptop_aucTure_list = [row.strip() for row in geptop_aucTure_file.readlines()]
	geptop_aucTure_dict = {}
	distance = {}
	distance_multi_aucTure = {}
	for row in spieces_distance_result:
		distance[row[0].split("/")[-1]] = row[1]
	for num in range(len(geptop_aucTure_list)):
		geptop_aucTure_dict["genomeDataset{}".format(num+1)] = float(geptop_aucTure_list[num])

	distance_sorted = sorted(distance.items(), key = lambda e:e[1],reverse=True)
	for num1 in range(5):
		spieces_name = distance_sorted[num1][0]
		distance_multi_aucTure["columns{}".format(num1+1)] = geptop_aucTure_dict[spieces_name]

	geotop_coefficient_data = pd.DataFrame(data=distance_multi_aucTure,index=["{}".format(userFileName)])
	geptop_coefficient = geptop_model.predict(geotop_coefficient_data)[0] if geptop_model.predict(geotop_coefficient_data)[0] <= 1 else 1.0
	
	scfm_coefficient_data = pd.DataFrame(data={"seed_5fold":seed_5folds_auc,"geptop_auc_pre":geptop_coefficient},index=["{}".format(userFileName)])
	scfm_coefficient = scfm_model.predict(scfm_coefficient_data)[0] if scfm_model.predict(scfm_coefficient_data)[0] <=1 else 1.0
   
	print("geptop_coefficient:{}".format(geptop_coefficient))
	print("scfm_coefficient:{}".format(scfm_coefficient))
	geptop_aucTure_file.close()

	return geptop_coefficient,scfm_coefficient

def result_final(geptop_id_score, scfm_id_score, geptop_coefficient,scfm_coefficient, Cutoff):
	"""Obtain the final results"""
	id_score_finally = {}

	for k,v in geptop_id_score.items():
		id_score_finally[k] = geptop_id_score[k]*(geptop_coefficient**26)+scfm_id_score[k]*(scfm_coefficient**26)
	score_finally = id_score_finally.values()
	max_score = max(score_finally)
	min_score = min(score_finally)
	EssentialGeneNumber = 0
	for k,v in id_score_finally.items():
		id_score_finally[k] = (id_score_finally[k]-min_score)/(max_score-min_score)
		if id_score_finally[k] >= Cutoff:
			EssentialGeneNumber += 1

	return id_score_finally, EssentialGeneNumber

def write_result(userFileName, Prediction, EssentialGeneNumber, Cutoff):
	"""Write the results into a file"""

	finalResultFile = open(rootdir+"/resultFile/{}.txt".format(userFileName.split('.')[0]), "w")
	finalResultFile.write("#Total %d genes are submitted\n" % len(Prediction))
	print("#Total %d genes are submitted\n" % len(Prediction))

	finalResultFile.write("#%d of them are predicted as essential genes\n" % EssentialGeneNumber)
	print("#%d of them are predicted as essential genes\n" % EssentialGeneNumber)

	finalResultFile.write('Class(essential gene:1,others:0)\tEssentiality Score\tID\n')
	print('Class(essential gene:1,others:0)\tEssentiality Score\tProtein\n')

	for k, v in Prediction.items():
		if v>Cutoff:
			finalResultFile.write('1\t\t\t\t\t\t\t'+str(float('%.4f'% v))+'\t\t\t\t\t'+k+'\n')
		else:
			finalResultFile.write('0\t\t\t\t\t\t\t'+str(float('%.4f'% v))+'\t\t\t\t\t'+k+'\n')
	finalResultFile.close()
	print("write successfulldy")


def workflow1(Target, userFileName):
	"""Workflow of homology_mapping"""

	blastdb = 'makeblastdb -in ' + Target+' -dbtype prot -parse_seqids -hash_index -out ' + os.path.join("{}/blastdatabase".format(rootdir),userFileName)
	print(blastdb)
	os.system(blastdb)

	Target_Record = SeqIO.parse(open(Target),'fasta')
	Geptop_prediction, GeneID, EssentialGeneid = information_exstract(Target_Record)

	TargetCV=CompositionVector(Target)

	p=Pool(8)  # Number of parallel processes
	result=[] 
	for root,dirs,files in os.walk('{}/datasets3'.format(rootdir)):
		for f in files:	 
			OGfaa=os.path.join(root,f)
			result.append(p.apply_async(homology_mapping,args=(TargetCV,Target,OGfaa,)))
	p.close()
	p.join()   

	spieces_distance_result = homo_result_extract(result, Geptop_prediction, EssentialGeneid,)
	Geptop_result = result1(Geptop_prediction, GeneID)
	remove_file(userFileName)
	print("workflow1 done")

	return spieces_distance_result, Geptop_result

def workflow2(Geptop_result, userFileName):
	"""Seed selection and extraction of sequence features"""
	seed_essential, seed_nonessential = seed_select(Geptop_result)
	id_seq = id_seq_get(userFileName)
	parallel_process = 8
	feature_result = multi_process_feature_get(id_seq, parallel_process)
	seed_feature = seed_feature_extract(feature_result,seed_essential, seed_nonessential)
	print("workflow2 done")

	return feature_result, seed_feature

def workflow3(feature, seed_feature, Geptop_result, spieces_distance_result, userFileName, Cutoff):
	""""""
	X = seed_feature.iloc[:,1:4044]
	Y = seed_feature['essensiality']
	smo = SMOTE(random_state = 608)
	X_ban, Y_ban = smo.fit_resample(X,Y)
	model_trained = train_model(X_ban, Y_ban)
	model_trained.fit(X_ban, Y_ban)
	seed_5folds_auc = cross_val_score(model_trained, X, Y, scoring='roc_auc', cv=5).mean()
	scfm_id_score = prediction(model_trained, feature)
	geptop_coefficient,scfm_coefficient = coefficient_caculate(userFileName, seed_5folds_auc, spieces_distance_result)
	result_finally, EssentialGeneNumber =  result_final(Geptop_result, scfm_id_score, geptop_coefficient,scfm_coefficient, Cutoff)
	write_result(userFileName, result_finally, EssentialGeneNumber, Cutoff)
	print("workflow3 done")


def main():
	""""""
	Target, Cutoff = getargu()
	userFileName = os.path.split(Target)[-1]
	translatedFile = translation(Target, userFileName)
	spieces_distance_result, Geptop_result, = workflow1(translatedFile, userFileName)
	feature, seed_feature = workflow2(Geptop_result, userFileName)
	workflow3(feature, seed_feature, Geptop_result, spieces_distance_result, userFileName, Cutoff)


if __name__ =='__main__':
	rootdir = "/****/***/Geptop3"  # Folder directory of Geptop3
	main()
