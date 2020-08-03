# Input/Sim {
import pandas as pd
# }

# Sim {
t_1 = ''
teams = ['Air Force', 'Akron', 'Alabama', 'Appalachian State', 'Arizona', 'Arizona State', 'Arkansas', 'Arkansas State', 'Army', 'Auburn', 'Ball State', 'Baylor','Boise State', 'Boston College', 'Bowling Green State', 'Brigham Young', 'Buffalo', 'California', 'Central Michigan', 'Charlotte', 'Cincinnati', 'Clemson', 'Coastal Carolina', 'Colorado', 'Colorado State', 'Connecticut', 'Duke', 'East Carolina', 'Eastern Michigan', 'Florida', 'Florida Atlantic', 'Florida International', 'Florida State', 'Fresno State', 'Georgia', 'Georgia Southern', 'Georgia State', 'Georgia Tech', 'Hawaii', 'Houston', 'Illinois', 'Indiana', 'Iowa', 'Iowa State', 'Kansas', 'Kansas State', 'Kent State', 'Kentucky', 'Louisiana', 'Louisiana Tech', 'Louisiana-Monroe', 'Louisville', 'LSU', 'Marshall', 'Maryland', 'Massachusetts', 'Memphis', 'Miami (FL)', 'Miami (OH)', 'Michigan', 'Michigan State', 'Middle Tennessee State', 'Minnesota', 'Mississippi State', 'Missouri', 'Navy', 'Nebraska', 'Nevada', 'Nevada-Las Vegas', 'New Mexico', 'New Mexico State', 'North Carolina', 'North Carolina State', 'North Texas', 'Northern Illinois', 'Northwestern', 'Notre Dame', 'Ohio', 'Ohio State', 'Oklahoma', 'Oklahoma State', 'Old Dominion', 'Ole Miss', 'Oregon', 'Oregon State', 'Penn State', 'Pitt', 'Purdue', 'Rice', 'Rutgers', 'San Diego State', 'San Jose State', 'SMU', 'South Carolina', 'South Florida', 'Southern Mississippi', 'Stanford', 'Syracuse', 'Temple', 'Tennessee', 'Texas', 'Texas A&M', 'Texas Christian', 'Texas State', 'Texas Tech', 'Toledo', 'Troy', 'Tulane', 'Tulsa', 'UAB', 'UCF', 'UCLA', 'USC', 'Utah', 'Utah State', 'UTEP', 'UTSA', 'Vanderbilt', 'Virginia', 'Virginia Tech', 'Wake Forest', 'Washington', 'Washington State', 'West Virginia', 'Western Kentucky', 'Western Michigan', 'Wisconsin', 'Wyoming']
p_score_r = 0
p_score_w = 0
p_score = 0
r_score_r = 0
r_score_w = 0
r_score = 0
f1_score = 0
cfb_cr = {}
for t_1 in teams:
    
    # Input/Sim {
    s = 13

    year = 2019
    y = str(year)
    yp = str(year - 1)
    # }

    # Sim {
    w = 0
    l = 0
    g = 1
    for g in range(g, s, 1):

        # Sim {
        t_1_m = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_1+y+'.csv',usecols=[7])
        t_nm = t_1_m.iat[g-1,0]
        if t_nm == 'Non-Major':
            continue

        t_1_s = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_1+y+'.csv',usecols=[6])
        t_2_l = t_1_s.iat[g-1,0]

        t_2_s = t_2_l.replace('(1)', '').replace('(2)', '').replace('(3)', '').replace('(4)', '').replace('(5)', '').replace('(6)', '').replace('(7)', '').replace('(8)', '').replace('(9)', '').replace('(10)', '').replace('(11)', '').replace('(12)', '').replace('(13)', '').replace('(14)', '').replace('(15)', '').replace('(16)', '').replace('(17)', '').replace('(18)', '').replace('(19)', '').replace('(20)', '').replace('(21)', '').replace('(22)', '').replace('(23)', '').replace('(24)', '').replace('(25)', '')
        if len(t_2_s) == len(t_2_l):
            t_2 = t_2_l
        else:
            t_2 = t_2_s[1:]

        h_1_s = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_1+y+'.csv',usecols=[5])
        h = h_1_s.iat[g-1,0]
        if h is '@':
            h_1 = 0
            h_2 = 3.5
        elif h is 'N':
            h_1 = 0
            h_2 = 0
        else:
            h_1 = 3.5
            h_2 = 0
        # }

        # Input/Sim {
        srsp = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/SRS'+yp+'.csv',usecols=[1, 8])
        srsp_1 = srsp.loc[srsp['School']==str(t_1), 'SRS']
        srsp_2 = srsp.loc[srsp['School']==str(t_2), 'SRS']
        srs = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/SRS'+y+'.csv',usecols=[1, 8])
        srs_1 = srs.loc[srsp['School']==str(t_1), 'SRS']
        srs_2 = srs.loc[srsp['School']==str(t_2), 'SRS']

        m_1 = (srsp_1 * ((s-g)/s) + srs_1 * (g/s))
        m_2 = (srsp_2 * ((s-g)/s) + srs_2 * (g/s))
        mm_sr = m_1.append(m_2)
        m = (mm_sr.iloc[0] - mm_sr.iloc[1]) / 2
        # }

        # Input/Sim {
        pyp_1 = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_1+yp+'.csv',usecols=[9])
        py_1 = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_1+y+'.csv',usecols=[9])
        oyp_1 = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_1+yp+'.csv',usecols=[10])
        oy_1 = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_1+y+'.csv',usecols=[10])
        pyp_2 = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_2+yp+'.csv',usecols=[9])
        py_2 = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_2+y+'.csv',usecols=[9])
        oyp_2 = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_2+yp+'.csv',usecols=[10])
        oy_2 = pd.read_csv('/Users/Dallon/Desktop/Python_Stuff/Sports_App/CFB_Team_Data/'+t_2+y+'.csv',usecols=[10])

        mmpyp_1 = pyp_1.mean() * ((s-g)/s) + pyp_1.median() * (g/s)
        sdpyp_1 = pyp_1.std()
        mmpy_1 = py_1.mean() * ((s-g)/s) + py_1.median() * (g/s)
        sdpy_1 = py_1.std()
        mmoyp_1 = oyp_1.mean() * ((s-g)/s) + oyp_1.median() * (g/s)
        sdoyp_1 = oyp_1.std()
        mmoy_1 = oy_1.mean() * ((s-g)/s) + oy_1.median() * (g/s)
        sdoy_1 = oy_1.std()
        mmpyp_2 = pyp_2.mean() * ((s-g)/s) + pyp_2.median() * (g/s)
        sdpyp_2 = pyp_2.std()
        mmpy_2 = py_2.mean() * ((s-g)/s) + py_2.median() * (g/s)
        sdpy_2 = py_2.std()
        mmoyp_2 = oyp_2.mean() * ((s-g)/s) + oyp_2.median() * (g/s)
        sdoyp_2 = oyp_2.std()
        mmoy_2 = oy_2.mean() * ((s-g)/s) + oy_2.median() * (g/s)
        sdoy_2 = oy_2.std()

        mmp_1 = mmpyp_1 * ((s-g)/s) + mmpy_1 * (g/s)
        sdp_1 = sdpyp_1 * ((s-g)/s) + sdpy_1 * (g/s)
        mmo_1 = mmoyp_1 * ((s-g)/s) + mmoy_1 * (g/s)
        sdo_1 = sdoyp_1 * ((s-g)/s) + sdoy_1 * (g/s)
        mmp_2 = mmpyp_2 * ((s-g)/s) + mmpy_2 * (g/s)
        sdp_2 = sdpyp_2 * ((s-g)/s) + sdpy_2 * (g/s)
        mmo_2 = mmoyp_2 * ((s-g)/s) + mmoy_2 * (g/s)
        sdo_2 = sdoyp_2 * ((s-g)/s) + sdoy_2 * (g/s)

        ms = mmp_1.append(mmo_2).append(sdp_1).append(sdo_2).append(mmp_2).append(mmo_1).append(sdp_2).append(sdo_1)
        # }

        # Sim {
        msh_1 = (ms.iloc[0] + ms.iloc[2] * 2) * .5 + (ms.iloc[1] + ms.iloc[3] * 2) * .5 + m + h_1
        msl_1 = (ms.iloc[0] - ms.iloc[2] * 2) * .5 + (ms.iloc[1] - ms.iloc[3] * 2) * .5 + m + h_1
        if msl_1 < 0:
            msl_1 = 0
        else:
            msl_1 = (ms.iloc[0] - ms.iloc[2] * 2) * .5 + (ms.iloc[1] - ms.iloc[3] * 2) * .5 + m + h_1
        msh_2 = (ms.iloc[4] + ms.iloc[6] * 2) * .5 + (ms.iloc[5] + ms.iloc[7] * 2) * .5 - m + h_2
        msl_2 = (ms.iloc[4] - ms.iloc[6] * 2) * .5 + (ms.iloc[5] - ms.iloc[7] * 2) * .5 - m + h_2
        if msl_2 < 0:
            msl_2 = 0
        else:
            msl_2 = (ms.iloc[4] - ms.iloc[6] * 2) * .5 + (ms.iloc[5] - ms.iloc[7] * 2) * .5 - m + h_2
        c = (msh_1 - msl_1 + 1) * (msh_2 - msl_2 + 1)
        # }

        # Input/Sim {
        n_1 = 0
        n_2 = 0
        if msh_2 > msh_1:
            n_2 = (msh_2 - msh_1) * (msh_1 - msl_1 + 1) + ((((msh_1 - msl_1 + 1) - 1) / 2) * (msh_1 - msl_1 + 1))
            n_1 = c - n_2 - msh_1
        elif msh_2 < msh_1:
            n_1 = (msh_1 - msh_2) * (msh_2 - msl_2 + 1) + ((((msh_2 - msl_2 + 1) - 1) / 2) * (msh_2 - msl_2 + 1))
            n_2 = c - n_1 - msh_2

        p_1 = (n_1 / c) * 100 / .95
        p_2 = (n_2 / c) * 100 / .95
        if p_1 > 99:
            p_1 = 99
            p_2 = 1
        elif p_2 > 99:
            p_2 = 99
            p_1 = 1

        m_1 = mmp_1.combine(mmo_2, max, fill_value=0)
        m_2 = mmo_1.combine(mmp_2, max, fill_value=0)
        # }

        # Sim {
        mm_1 = m_1.mean() + m + h_1
        mm_2 = m_2.mean() - m + h_2
        #mm = - mm_1 + mm_2
        mm = - p_1 + p_2

        po = - py_1.iat[g-1, 0] + oy_1.iat[g-1, 0]
        moe =  mm - po

        if po < 0 and mm < 0:
            p_score_r += 1
        elif po > 0 and mm > 0:
            r_score_r += 1
        elif po < 0 and mm > 0:
            r_score_w += 1
        elif po > 0 and mm < 0:
            p_score_w += 1


        if mm < 0:
            w += 1
        elif mm > 0:
            l += 1

        print(g, t_1, int(p_1), round(mm, 1), po, t_2, int(p_2), sep=',')
        g += 1
        # }

    cfb_cr[str(t_1)] = (int(w), int(l))
    # }

print(cfb_cr)
p_score = (p_score_r / (p_score_r + p_score_w))
r_score = (r_score_r / (r_score_r + r_score_w))
#f1_score = p_score * .5 + r_score * .5
f1_score = ((p_score_r + r_score_r) / (p_score_r + p_score_w + r_score_r + r_score_w))
print('f1_score', round(f1_score, 4), sep=',')
# }