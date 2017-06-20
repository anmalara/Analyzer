import FWCore.ParameterSet.Config as cms

import sys
print 1

files=['root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/007A655B-7FB7-E611-A316-0025905B858C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/02392142-80B7-E611-8A2F-0CC47A7C351E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/0896780D-7DB7-E611-9597-002618FDA204.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/0A5A5DE8-83B7-E611-A144-0CC47A4D76CC.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/0CC0527C-7EB7-E611-846A-0CC47A78A446.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/0EE985DF-7BB7-E611-880E-0CC47A78A3D8.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/1821BBCA-82B7-E611-895B-0CC47A7C351E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/1A09A4C8-7FB7-E611-BD49-0025905B859A.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/1A602152-81B7-E611-8A13-0CC47A7AB7A0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/1A98CDE8-83B7-E611-848B-0025905A6080.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/20537522-82B7-E611-ABA8-0CC47A4D767E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/22EBE3DC-7DB7-E611-B070-0CC47A4C8E5E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/24B31858-81B7-E611-BFC0-0025905B8568.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/2E3318AA-82B7-E611-9E31-0CC47A4C8E70.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/36C6EBCC-82B7-E611-8C94-0CC47A7AB7A0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/3E0A813B-83B7-E611-A3D6-0CC47A7C35B2.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/3EF0E54A-88B7-E611-8778-0CC47A4D76CC.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/4861B1DE-7BB7-E611-B15F-0CC47A4C8E70.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/4A1ACDB1-7CB7-E611-AE3E-0025905B861C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/4AB933C3-7EB7-E611-9C2F-0CC47A4D76CC.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/4CB87587-8BB7-E611-A81E-0CC47A78A41C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/50C56517-7DB7-E611-8AB1-0025905B85D0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/50ED67A1-7DB7-E611-A998-0CC47A4C8E70.root',
####'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/56317534-84B7-E611-9FD5-0CC47A7C351E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/5AF4294C-7DB7-E611-BF92-0CC47A78A41C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/5E57D496-80B7-E611-A953-0CC47A4C8E5E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/60EAAA21-7BB7-E611-876C-0CC47A78A478.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/62A6E890-83B7-E611-B6A6-0CC47A4C8E5E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/62CA8098-83B7-E611-B623-0CC47A78A446.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/62DBCCA4-83B7-E611-B2A9-0CC47A4C8F0C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/647BCF80-80B7-E611-ABE2-0CC47A4D764C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/689F5730-85B7-E611-89F0-0CC47A7C35B2.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/6A186EC2-7FB7-E611-AECE-0CC47A4D7630.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/6AAB74AE-7CB7-E611-89F5-0CC47A7C3420.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/6AC344E0-80B7-E611-9BFE-0CC47A7C35F4.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/6ED0A541-80B7-E611-8238-0CC47A78A446.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/722073BA-7CB7-E611-B0DD-0CC47A4C8E7E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/7E5AE41F-83B7-E611-A9C7-0025905B860C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/88503387-8BB7-E611-BE79-0025905A6064.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/8E3EE9E0-82B7-E611-AC03-0025905A60FE.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/92E4F73E-80B7-E611-8E2A-0CC47A4C8F0C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/94933B97-7EB7-E611-90E4-0025905B85D0.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/9E26BFFF-81B7-E611-A7EA-0025905B85BE.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/A4D81D22-7BB7-E611-8375-0CC47A4C8E70.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/A6DFCDAD-7CB7-E611-ACD8-0CC47A4C8E98.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/A81D6158-81B7-E611-BD83-0025905B8566.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/AEA8BC37-83B7-E611-A5B3-0CC47A78A4BA.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/B22C4221-82B7-E611-BDF8-0025905A60CA.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/B2C0D123-7FB7-E611-8AF5-0025905B855E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/B64C8853-81B7-E611-83E3-0CC47A7C3420.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/B6A49859-81B7-E611-90A4-0025905B860C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/BA5F6C11-7DB7-E611-99F8-0CC47A7C35B2.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/C03C47E4-7DB7-E611-AA3D-0CC47A4D764C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/C4588417-7FB7-E611-81D9-0025905B861C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/C4C65B0B-7DB7-E611-8910-0CC47A4D7602.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/C89CBCEA-80B7-E611-8A6B-0CC47A4C8E70.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/CADD3AC2-7EB7-E611-ABA7-0CC47A7C35B2.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/CEBDEA4B-84B7-E611-96BF-0025905B85BE.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/CEC9B44D-7FB7-E611-9127-0CC47A7C35F4.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/D2809F1F-85B7-E611-8896-0025905A60FE.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/D2DD0252-83B7-E611-9C5E-0025905B8566.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/DA8BB8BA-7CB7-E611-8DDA-0025905A608A.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/DA944621-7FB7-E611-BA44-0025905B8566.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/DAB9E887-83B7-E611-BB2C-0025905B855E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/DC62889D-7DB7-E611-B8C8-0025905A6080.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/E02E3FF9-81B7-E611-9C19-0CC47A4D764C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/E4A9A373-7EB7-E611-9283-0CC47A4C8F0C.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/E4BAF56F-83B7-E611-ADAD-0025905B859A.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/EAF2C19B-7FB7-E611-BA51-0025905A6080.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/F04AE967-88B7-E611-A822-0CC47A4D7630.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/F20E0EF7-81B7-E611-9E84-0CC47A4C8E5E.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/F2E86A1F-82B7-E611-8A44-0025905A6080.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/F419E2DF-7BB7-E611-8EC4-0CC47A4C8E82.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/F4553951-81B7-E611-9F8B-0CC47A7452DA.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/FAC73FFC-89B7-E611-BD93-0025905B859A.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/FAE96376-83B7-E611-8DAB-0CC47A4D7674.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/FC38391F-82B7-E611-99A2-0025905A6080.root',
'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/FE106C51-7FB7-E611-8A02-0025905B85BE.root']

print len(files)
print files[int(sys.argv[2]):(int(sys.argv[3]))]


process = cms.Process("Pippo")

#process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_38T_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_v3'


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        files[int(sys.argv[2]):(int(sys.argv[3]))]
        #'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/007A655B-7FB7-E611-A316-0025905B858C.root'
        #'file:file_QCD_Pt470to600.root'
        #'file:file_QCD_Pt120to170.root'
        #'file:file_QCD_Pt1800to2400.root'
      #'file:/scratch/lgiannini/FAE96376-83B7-E611-8DAB-0CC47A4D7674.root',




        
        #'root://xrootd-cms.infn.it//store/relval/CMSSW_8_1_0_pre16/RelValTTbar_13/GEN-SIM-RECO/PU25ns_81X_upgrade2017_realistic_v22_HS-v1/10000/007A655B-7FB7-E611-A316-0025905B858C.root'
        #'root://xrootd-cms.infn.it//store/mc/RunIISpring16MiniAODv2/TT_TuneCUETP8M1_13TeV-powheg-pythia8-evtgen/MINIAODSIM/PUSpring16_80X_mcRun2_asymptotic_2016_miniAODv2_v0-v1/00000/02EA5736-2C3F-E611-AAD2-0CC47A4D99A4.root'
    
    )#,
    #skipEvents=cms.untracked.uint32(36590)
)

#for MC jet flavour
process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
process.load("PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi")
process.ak4JetFlavourInfos.jets = cms.InputTag("ak4PFJetsCHS")
process.ak4JetFlavourInfos.hadronFlavourHasPriority = cms.bool(True)
process.flavourSeq = cms.Sequence(
    process.selectedHadronsAndPartons *
    process.ak4JetFlavourInfos
)


process.analyzer1 = cms.EDAnalyzer('AnalyzerSignedIP_RECO_wArr',
 jetMCSrc = cms.InputTag("ak4JetFlavourInfos"),
#    bits = cms.InputTag("TriggerResults","","HLT"),
#    prescales = cms.InputTag("patTrigger"),
#    objects = cms.InputTag("selectedPatTrigger"),
#    generatorInfo = cms.InputTag("generator"),
)

outfile='recoFiles/recoArrays_dd'+sys.argv[2]+'.root'

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(outfile) )


process.p = cms.Path(process.flavourSeq+process.analyzer1)


