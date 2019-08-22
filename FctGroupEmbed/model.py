# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 13:03:16 2019

@author: simoussi
"""


from keras.models import Model
from keras.layers import Input, Dense
import pandas as pd
import tensorflow as tf
from sklearn.model_selection import train_test_split
import os
from keras import backend as KE

def r2_keras(y_true, y_pred):
    SS_res =  KE.sum(KE.square(y_true - y_pred)) 
    SS_tot = KE.sum(KE.square(y_true - KE.mean(y_true))) 
    return ( 1 - SS_res/(SS_tot + KE.epsilon()) )

class FctGroupModel():
    def __init__(self,archi_config,learn_confuig,pool):
        self.archi_config=archi_config
        self.learn_config=learn_config
        self.names=pool
        
    def create_model(self):
        print("Creating model architecture")
        input_layer=Input(shape=(archi_config.get('inSize'),),dtype=tf.float32,name="taxo_composition")
        hidden=Dense(archi_config.get('embSize'),activation=archi_config.get('hidAct'),kernel_regularizer=archi_config.get('reg'))(input_layer)
        output=Dense(archi_config.get('outSize'))(hidden)
        
        self.model=Model(inputs=input_layer,outputs=output)
        
    def compile_model(self):
        self.model.compile(optimizer=self.learn_config.get('optim'),loss=self.learn_config.get('loss'),metrics=self.learn_config.get('metric'))
    
    def train_model(self,dataset):
        history=self.model.fit(x=dataset.get('X'),y=dataset.get('Y'),batch_size=self.learn_config.get('batchSize'),epochs=self.learn_config.get('epoch'),validation_split=0.1,verbose=1)
        return history
    
    def evaluate_model(self,dataset):
        eval_test=self.model.evaluate(x=dataset.get('X'),y=dataset.get('Y'))
        return eval_test
    
    def save_embeddings(self,save=True,file="embeddings.tsv"):
        W=pd.DataFrame(data=self.model.get_weights()[0])
        W.to_csv(file+".tsv",sep="\t",index=False,decimal=".",header=False)
        
        pd.DataFrame(self.names).to_csv(file+"_metadata.tsv",sep="\t",index=False,decimal=".",header=False)
 
    
def main(data_file,out_cols,archi_config,learn_config,save=True,emb_file="emb.tsv"):
    data=pd.read_csv(data_file,sep=",",decimal=".")
    
    fct=data[out_cols]
    composition=data.drop(out_cols,axis=1)
    pool=composition.columns.tolist()
    
    n,p=composition.shape
    m=fct.shape[1]
    
    archi_config.update(dict(inSize=p,outSize=m))
    fgm=FctGroupModel(archi_config,learn_config,pool)
    fgm.create_model()
    fgm.model.summary()
    fgm.compile_model()
    
    X_train, X_test, y_train, y_test = train_test_split(
     composition, fct, test_size=0.2, random_state=42)
    
    #### Begin training ###
    train_dataset=dict(X=X_train.values.astype('float'),Y=y_train.iloc[:,0:m].values)
    test_dataset=dict(X=X_test.values.astype('float'),Y=y_test.iloc[:,0:m].values)
    
    hist=fgm.train_model(train_dataset)
    del hist
    score=fgm.evaluate_model(test_dataset)
    print("Score = ",score)
    
    if save:
        fgm.save_embeddings(save,emb_file)
    
    return fgm, score

if __name__== "__main__":
    data_file="data/tilman.csv"
    scores=[]
    for K in [8]:#[2,4,8]:
        print("Training for K = ",K)
        out_dir="output/"+str(K)+"/"
        os.makedirs(out_dir,exist_ok=True)
        
        for opt in ["sgd"]:#["adam","rmsp","sgd","adagrad"]:
            for act in ["relu"]:#["relu",None]:
                archi_config={
                    "embSize":K,
                    "reg":"l2",
                    "outSize":1,
                    "hidAct":act,
                    }
                
                learn_config={
                        "batchSize":8,
                        "epoch":50,
                        "optim":"adam",
                        "loss":"MSE",
                        "metric":["MAE",r2_keras]
                        }
                
                if act is None:    
                    fgm, score=main(data_file,['Fobs'],archi_config,learn_config,save=True,emb_file=out_dir+opt+"_none_embeddings")
                else:
                    fgm, score=main(data_file,['Fobs'],archi_config,learn_config,save=True,emb_file=out_dir+opt+"_"+act+"_embeddings")
                
                scores.append({
                        'K':K,
                        'opt':opt,
                        'act':act,
                        'mae_score':score[1],
                        'mse_score':score[0]
                        })
                
    
    df=pd.DataFrame.from_dict(scores)
        
    
    
