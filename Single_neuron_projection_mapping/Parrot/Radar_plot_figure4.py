import numpy as np
import matplotlib.pyplot as plt

# 1. Define brain area
# 
#labels = ['NAO', 'MO', 'LPOm', 'DMm','DM', 'nXIIts', 'PAM','NAOc', 'MOc', 'LPOmc', 'DMmc','DMc', 'nXIItsc', 'PAMc']

labels = ['NAO', 'MO', 'LPOm', 'DMm','DM', 'nXIIts', 'PAM', 'PAMc', 'nXIItsc' ,'DMc', 'DMmc', 'LPOmc', 'MOc', 'NAOc']
num_vars = len(labels)



# normalize to nXIIts
zebra_finch_proportions = [0, 0, 0, 4.3520714, 28.5953411, 39.2662854, 27.7863021, 0, 0, 0, 0, 0, 0, 0]
#parrot_proportions = [48.16414687, 51.90784737, 20.08639309, 15.55075594, 14.11087113, 43.26853852, 27.06983441, 43.26853852, 55.36357091, 24.550036, 4.391648668, 6.623470122, 42.40460763, 20.87832973]

parrot_proportions = [48.16414687, 51.90784737, 20.08639309, 15.55075594, 14.11087113, 43.26853852, 27.06983441, 20.87832973, 42.40460763, 6.623470122, 4.391648668, 24.550036, 55.36357091,  43.26853852]


def close_loop(data):
    return np.concatenate((data, [data[0]]))

zebra_finch_plot_data = close_loop(zebra_finch_proportions)
parrot_plot_data = close_loop(parrot_proportions)


angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1] 

# 2. draw the picture
fig, ax = plt.subplots(figsize=(9, 9), subplot_kw=dict(polar=True))


ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)


ax.plot(angles, zebra_finch_plot_data, color='#2E3192', linewidth=2.5, label='Zebra Finch(RA)', marker='o', markersize=4)
ax.fill(angles, zebra_finch_plot_data, color='#2E3192', alpha=0.2)


ax.plot(angles, parrot_plot_data, color='#C1272D', linewidth=2.5, label='Parrot (AAcd)', marker='o', markersize=4)
ax.fill(angles, parrot_plot_data, color='#C1272D', alpha=0.2)


ax.set_thetagrids(np.degrees(angles[:-1]), labels, fontsize=11, fontweight='bold')


max_val = 50 
ax.set_ylim(0, max_val)


ticks = [20, 40, 60]
ax.set_rticks(ticks)
ax.set_yticklabels([f"{t}%" for t in ticks], color="grey", size=9)
ax.set_rlabel_position(180)


ax.spines['polar'].set_visible(False)


plt.title('Brain Projection Profile Comparison\n(Percentage of Total Axonal Output)', 
          size=14, fontweight='bold', pad=40)

plt.legend(loc='upper right', bbox_to_anchor=(1.15, 1.1), fontsize=10, frameon=True)

plt.tight_layout()

# 3. Save picture
plt.savefig('parrot_zebra_brain_radar_chart_final.eps', dpi=300, bbox_inches='tight')
plt.show()