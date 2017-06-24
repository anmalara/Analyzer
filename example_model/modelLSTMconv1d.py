#from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Convolution2D, merge, Convolution1D, Reshape, LSTM, Permute
from keras.models import Model, load_model
from keras.layers import Input 
from random import randint
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

import numpy
seed = 7 
numpy.random.seed(seed)
import keras 

losses = []
accuracy = []
val_losses = []
val_accuracy = []

print "corrections"

f='../../DataMiniAODNew/jvars_420000_430000_MiniArrays0.npy'
jvars=numpy.load(f)
svars=numpy.load(f.replace('jvars_', 'svars'))
tvars=numpy.load(f.replace('jvars', 'tvars'))
jcol_mean = numpy.nanmean(jvars,axis=0)
jcol_mean[4]=0
jcol_std = numpy.nanstd(jvars,axis=0)
jcol_std[4]=1
scol_mean = numpy.nanmean(svars,axis=0)            
scol_std = numpy.nanstd(svars,axis=0)
tcol_mean = numpy.nanmean(tvars,axis=0)            
tcol_std = numpy.nanstd(tvars,axis=0)
jcol_mean=numpy.append(jcol_mean,[0,0,0,0])
jcol_std=numpy.append(jcol_std,[1,1,1,1])
#if normed==True:
#    col_mean = numpy.nanmean(a,axis=0)            
#    col_std = numpy.std(a,axis=0)
#    a[:]=a[:]-col_mean
#    a[:]=a[:]-col_std

def make_vars(f):
    jvars=numpy.load(f)
    svars=numpy.load(f.replace('jvarsNew_', 'svars'))
    tvars=numpy.load(f.replace('jvarsNew', 'tvars'))
    jinds = numpy.where(numpy.isnan(jvars))
    jevents_to_remove=numpy.unique(jinds[0])
    sinds = numpy.where(numpy.isnan(svars))
    sevents_to_remove=numpy.unique(sinds[0])
    tinds = numpy.where(numpy.isnan(tvars))
    tevents_to_remove=numpy.unique(tinds[0])

    events_to_remove=numpy.unique(numpy.concatenate((tevents_to_remove, sevents_to_remove, jevents_to_remove)))
    events_to_remove=numpy.sort(events_to_remove)
    events_to_remove=events_to_remove[::-1]

    print events_to_remove
    for ev in events_to_remove:            
        jvars=numpy.delete(jvars, ev,0)
        svars=numpy.delete(svars, ev,0)
        tvars=numpy.delete(tvars, ev,0)
        print "removed event:  ", ev

    print "mean and std corrections"
    print jvars[0] , svars[0][0]
    jvars[:]=jvars[:]-jcol_mean
    jvars[:]=jvars[:]/jcol_std
    svars[:]=svars[:]-scol_mean*(svars[:]!=0)
    svars[:]=svars[:]*(svars[:]!=0)/scol_std
    tvars[:]=tvars[:]-tcol_mean*(tvars[:]!=0)
    tvars[:]=tvars[:]*(tvars[:]!=0)/tcol_std
    print jvars[0],  svars[0][0]
    svars=numpy.swapaxes(svars, 1,2) 
    tvars=numpy.swapaxes(tvars, 1,2)
    tvars=tvars[:,:,19:739]
    tvars=tvars.reshape((len(tvars),10,20,36))
    #tvars=numpy.swapaxes(tvars, 1,3)
    return jvars, svars, tvars




remove_events=True
normed=False

print "inputs"
#Inputs=[Input(shape=(4,)) ]
#Inputs+=[Input(shape=(10,19))]
#Inputs+=[Input(shape=(20,36)) for i in range(10) ]


#dropoutRate=0.1

#tracks=Inputs[1]

#tracks  = Convolution1D(64, 1, init='lecun_uniform',  activation='relu')(tracks)
#tracks = Dropout(dropoutRate)(tracks)
#tracks  = Convolution1D(32, 1, init='lecun_uniform',  activation='relu')(tracks)
#tracks = Dropout(dropoutRate)(tracks)
#tracks  = Convolution1D(32, 1, init='lecun_uniform',  activation='relu')(tracks)
#tracks = Dropout(dropoutRate)(tracks)
#tracks  = Convolution1D(8, 1, init='lecun_uniform',  activation='relu')(tracks)
        
#tracks = Flatten()(tracks)

#seedLSTM=LSTM(64)

#seed0=seedLSTM(Inputs[2])
#seed1=seedLSTM(Inputs[3])
#seed2=seedLSTM(Inputs[4])
#seed3=seedLSTM(Inputs[5])
#seed4=seedLSTM(Inputs[6])
#seed5=seedLSTM(Inputs[7])
#seed6=seedLSTM(Inputs[8])
#seed7=seedLSTM(Inputs[9])
#seed8=seedLSTM(Inputs[10])
#seed9=seedLSTM(Inputs[11])

#v_seq = merge( [seed0, seed1, seed2, seed3, seed4, seed5, seed6, seed7, seed8, seed9] , mode='concat')
#v_seq = Reshape( (64,10))(v_seq)
#v_seq = Permute( (2,1))(v_seq)
#v_seq  = Convolution1D(64, 1, init='lecun_uniform',  activation='relu')(v_seq)
#v_seq = Dropout(dropoutRate)(v_seq)
#v_seq  = Convolution1D(32, 1, init='lecun_uniform',  activation='relu')(v_seq)
#v_seq = Dropout(dropoutRate)(v_seq)
#v_seq  = Convolution1D(32, 1, init='lecun_uniform',  activation='relu')(v_seq)
#v_seq = Dropout(dropoutRate)(v_seq)
#v_seq  = Convolution1D(8, 1, init='lecun_uniform',  activation='relu')(v_seq)

#v_seq = Flatten()(v_seq)

#x = merge( [Inputs[0] , tracks, v_seq] , mode='concat')


#x=  Dense(200, activation='relu',init='lecun_uniform')(x)
#x = Dropout(dropoutRate)(x)
#x=  Dense(100, activation='relu',init='lecun_uniform')(x)
#x = Dropout(dropoutRate)(x)
#x=  Dense(100, activation='relu',init='lecun_uniform')(x)
#x = Dropout(dropoutRate)(x)
#x=  Dense(100, activation='relu',init='lecun_uniform')(x)
#x = Dropout(dropoutRate)(x)
#x=  Dense(100, activation='relu',init='lecun_uniform')(x)
#x = Dropout(dropoutRate)(x)
#x=  Dense(100, activation='relu',init='lecun_uniform')(x)
#x = Dropout(dropoutRate)(x)
#predictions = Dense(1, activation='sigmoid',init='lecun_uniform')(x)
#model = Model(input=Inputs, output=predictions)
model=load_model("../LSTMconv1d_NewMCTruth/eval3Test_at_epoch9.h5")

print'compiling'

adam=keras.optimizers.Adam(lr=0.0001)
model.compile(loss='binary_crossentropy', optimizer=adam,metrics=['accuracy', 'fmeasure'])

model.summary()

import glob

files=glob.glob('../../DataMiniAODNew/*jvarsNew*')
print len(files)

val_files=glob.glob('../../DataMiniAODNewValidation/*jvarsNew*')
print len(val_files)
#val_eval=['../val_data/jvars_400000_410000_testMini_Arrays2.npy', '../val_data/jvars_410000_420000_testMini_Arrays2.npy', '../val_data/jvars_420000_430000_testMini_Arrays2.npy', '../val_data/jvars_390000_400000_testMini_Arrays2.npy',
#'../val_data/jvars_550000_560000_testMini_Arrays3.npy','../val_data/jvars_560000_570000_testMini_Arrays3.npy','../val_data/jvars_410000_420000_testMini_Arrays2.npy',  '../val_data/jvars_540000_550000_testMini_Arrays3.npy',
 #'../val_data/jvars_570000_580000_testMini_Arrays3.npy']


for epoch in range(10,100):

    for f in files: 
        
        jvars, svars, tvars=make_vars(f)
        
        print f

        data=[jvars[:,0:4]]
        data=data+[svars]
        data=data+[tvars[:,i,:,:] for i in range(10)]
        
        
        v_file=val_files[randint(0,50)]

        v_jvars, v_svars, v_tvars=make_vars(v_file )
#        v_jvars=numpy.load(v_file)
#        v_svars=numpy.load(v_file.replace('jvars_', 'svars'))
#        v_tvars=numpy.load(v_file.replace('jvars', 'tvars'))
#        
#        if remove_events==True:
#            del_nan(v_jvars, v_svars, v_tvars)
        
        print v_file    
        v_data=[v_jvars[:,0:4]]
        v_data=v_data+[v_svars]
        v_data=v_data+[v_tvars[:,i,:,:] for i in range(10)]
        
        print len(v_data), len(data)
        print "fitting"

        

        history=model.fit(data,jvars[:,4:5]==5, validation_data=(v_data,v_jvars[:,4:5]==5) , batch_size=512,  nb_epoch=1, verbose=1)
        
        losses.append(history.history['loss'][0])
        accuracy.append(history.history['acc'][0])
        val_losses.append(history.history['val_loss'][0])
        val_accuracy.append(history.history['val_acc'][0])
     
    plt.cla() 
    for v_file in val_files:
        v_jvars, v_svars, v_tvars=make_vars(v_file )
        v_data=[v_jvars[:,0:4]]
        v_data=v_data+[v_svars]
        v_data=v_data+[v_tvars[:,i,:,:] for i in range(10)]
        pred_Part=model.predict(v_data)
       
        try:
            v_jvars_all=numpy.concatenate((v_jvars_all,v_jvars),axis=0)
            predictions=numpy.concatenate((predictions,pred_Part),axis=0)
            print "try"
        except:
            print "except"
            v_jvars_all=v_jvars
            predictions=pred_Part
     
    print "v_Data", v_jvars_all.shape    
        
    out_file = open("log_roc_curves.txt","a")
    
    out_file.write("at epoch "+str(epoch)+"\n")
    print "at epoch", epoch
    
    fpr, tpr, _ = roc_curve(v_jvars_all[:,4:5]==5, predictions)    
    print "ROC AUC:  ",auc(fpr, tpr)
    out_file.write("ROC AUC:  "+str(auc(fpr, tpr))+"\n")
    print "no charmed hadrons"
    out_file.write("no charmed hadrons ")
    fpr, tpr, _ = roc_curve(v_jvars_all[v_jvars_all[:,4]!=4][:,4:5]==5, predictions[v_jvars_all[:,4]!=4])    
    print "ROC AUC:  ",auc(fpr, tpr)
    out_file.write("ROC AUC:  "+str(auc(fpr, tpr))+"\n")
    print "B vs D"
    out_file.write("B vs D ")        
    fpr, tpr, _ = roc_curve(v_jvars_all[(v_jvars_all[:,4]>=4)*(v_jvars_all[:,4]<=5)][:,4:5]==5, predictions[(v_jvars_all[:,4]>=4)*(v_jvars_all[:,4]<=5)]) 
    print "ROC AUC:  ",auc(fpr, tpr)
    out_file.write("ROC AUC:  "+str(auc(fpr, tpr))+"\n")
    print "B_DUSG"
    out_file.write("B_DUSG ")
    fpr, tpr, _ = roc_curve(v_jvars_all[(v_jvars_all[:,4]!=4)*(v_jvars_all[:,4]!=0)][:,4:5]==5, predictions[(v_jvars_all[:,4]!=4)*(v_jvars_all[:,4]!=0)])    
    print "ROC AUC:  ",auc(fpr, tpr)
    out_file.write("ROC AUC:  "+str(auc(fpr, tpr))+"\n")
    plt.semilogy(tpr, fpr, label='ROC curve')
    
    
    fpr, tpr, _ = roc_curve(v_jvars_all[(v_jvars_all[:,4]!=4)*(v_jvars_all[:,4]!=0)][:,4:5]==5, v_jvars_all[(v_jvars_all[:,4]!=4)*(v_jvars_all[:,4]!=0)][:,5])    
    print "ROC AUC:  ",auc(fpr, tpr)
    out_file.write("ROC AUC:  "+str(auc(fpr, tpr))+"\n")
    plt.semilogy(tpr, fpr, label='ROC curve')

    fpr, tpr, _ = roc_curve(v_jvars_all[(v_jvars_all[:,4]!=4)*(v_jvars_all[:,4]!=0)][:,4:5]==5, v_jvars_all[(v_jvars_all[:,4]!=4)*(v_jvars_all[:,4]!=0)][:,6])    
    print "ROC AUC:  ",auc(fpr, tpr)
    out_file.write("ROC AUC:  "+str(auc(fpr, tpr))+"\n")
    plt.semilogy(tpr, fpr, label='ROC curve')
    
        
    plt.semilogy([0, 1], [0, 1], 'k--')
    plt.xticks( [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1] )
    plt.grid(True, which='both')
    plt.xlim([0.0, 1.05])
    plt.ylim([0.0, 1.05])
    plt.xlabel('b-jet efficiency')
    plt.ylabel('b-jet mistag rate')
    plt.title('Evaluation ROC') 

    plt.savefig("eval3Test_at_epoch"+str(epoch)+".png")

    plt.cla()

       
    plt.cla()
    plt.plot(losses)
    plt.plot(accuracy)
    plt.plot(val_losses)
    plt.plot(val_accuracy)
    plt.savefig("eval3LossHistory_at_epoch"+str(epoch)+".png")

    model.save("eval3Test_at_epoch"+str(epoch)+".h5") 
    del v_jvars_all, predictions
    out_file.close()
 
