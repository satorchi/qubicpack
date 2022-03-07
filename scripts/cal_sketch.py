'''
$Id: cal_sketch.py
$auth: Steve Torchinsky <satorchi@apc.in2p3.fr>
$created: Fri 27 Apr 2018 07:30:58 CEST
$license: GPLv3 or later, see https://www.gnu.org/licenses/gpl-3.0.txt

          This is free software: you are free to change and
          redistribute it.  There is NO WARRANTY, to the extent
          permitted by law.

draw the optical path of the QUBIC calibration source at APC
'''

from math import *
import matplotlib.pyplot as plt
from qubicpack.utilities import figure_window_title
plt.close('all')
def cal_sketch(pointing_ang=50,     # QUBIC nominal pointing elevation angle
               D=9.0,               # horizontal distance to QUBIC from the calibration source
               Hqubic=2.571,        # vertical height of QUBIC window from the floor
               Hmirror=3.7,         # vertical height of the centre of the flat mirror
               Hsource=4.0,         # vertical height of the calibration source
               mirror_len=1.4,      # length of the mirror
               cryo_width=1.436752, # width of the QUBIC cryostat
               cryo_height=1.542618,# height of the QUBIC cryostat
               win_offset=0.241979, # distance of the window from the centre axis of the cryostat
               xlength=13,          # x extent of the drawing
               ylength=7,           # y extent of the drawing
               transparency=True,
               figsize=(20,16),
               fontsize=16):
    
    window_pos=[-D,Hqubic]
    cal_pos=[0,Hsource]

    pointing_ang_rad=radians(pointing_ang)
    
    fig=plt.figure(figsize=figsize)
    figure_window_title(fig,'Calibration Source Sketch')
    ax = fig.add_axes((0.05,0.05,0.90,0.90))
    


    # QUBIC window marker
    x=window_pos[0]
    y=window_pos[1]
    ax.plot(x,y,color='red',linestyle='none',marker='D')
    txt='   QUBIC position (%.2fm, %.2fm)' % (D,Hqubic)
    ax.text(x+0.7,y,txt,color='black',fontsize=fontsize,ha='left',va='bottom')

    # cal source position marker
    x=cal_pos[0]
    y=cal_pos[1]
    ax.plot(x,y,color='green',linestyle='none',marker='D')
    txt='cal source height %.2fm  ' % cal_pos[1]
    ax.text(x,y,txt,color='black',fontsize=fontsize,ha='right',va='bottom')

    # line showing QUBIC line of sight
    s=2
    hyp=s/sin(pointing_ang_rad)
    c=hyp*cos(pointing_ang_rad)
    pt1=window_pos
    pt2=[window_pos[0]-c,window_pos[1]+s]

    # intersection of cal mirror height with QUBIC line of sight
    # y = mx+b : line of sight
    # y = Hmirror
    m=(pt2[1]-pt1[1])/(pt2[0]-pt1[0])
    b=pt1[1] - m*pt1[0]

    xmirror=(Hmirror-b)/m
    mirror_pos=[xmirror,Hmirror]

    x=[window_pos[0],mirror_pos[0]]
    y=[window_pos[1],mirror_pos[1]]
    ax.plot(x,y,color='red',linewidth=3)

    # line from cal source to mirror
    x=[cal_pos[0],mirror_pos[0]]
    y=[cal_pos[1],mirror_pos[1]]
    ax.plot(x,y,color='green',linewidth=3)

    delta_s=y[0]-y[1]
    delta_c=-(x[1]-x[0])
    delta_ang_rad=atan(delta_s/delta_c)
    delta_ang=degrees(delta_ang_rad)

    # reflection angle at the mirror
    reflection_ang=delta_ang+pointing_ang

    # mirror has normal at half the reflection angle from the QUBIC pointing angle
    mirror_norm_ang=pointing_ang-0.5*reflection_ang
    mirror_norm_ang_rad=radians(mirror_norm_ang)
    mirror_ang=mirror_norm_ang+90
    mirror_ang_rad=radians(mirror_ang)

    # a dashed line to show the mirror normal
    L=1.0
    dx=L*cos(mirror_norm_ang_rad)
    dy=L*sin(mirror_norm_ang_rad)
    x2=mirror_pos[0]+dx
    y2=mirror_pos[1]-dy
    x=[mirror_pos[0],x2]
    y=[mirror_pos[1],y2]
    ax.plot(x,y,color='black',linestyle='dashed',linewidth=1)
    txt='reflection angle at mirror = %.2f$^\circ$' % reflection_ang
    ax.text(x2,y2,txt,color='black',fontsize=fontsize,ha='left',va='center')
    
    # a line 1.4m long to show the mirror
    mirror_dx=0.5*mirror_len*cos(mirror_ang_rad)
    mirror_dy=0.5*mirror_len*sin(mirror_ang_rad)

    x=[mirror_pos[0]+mirror_dx,mirror_pos[0]-mirror_dx]
    y=[mirror_pos[1]-mirror_dy,mirror_pos[1]+mirror_dy]
    ax.plot(x,y,color='black',linewidth=3)

    # a dashed line showing the mirror horizontal extension
    x=[mirror_pos[0]-mirror_dx,mirror_pos[0]+mirror_dx]
    y=[mirror_pos[1]+mirror_dy,mirror_pos[1]+mirror_dy]
    ax.plot(x,y,color='black',linestyle='dashed')
    #ax.plot([x[1],x[1]],[y[1],0],color='black',linestyle='dashed')

    txt='mirror horizontal extension = %.2fm' % abs(2*mirror_dx)
    txt='%.2fm' % abs(2*mirror_dx)
    ax.text(x[0],y[0],txt,ha='right',va='bottom',fontsize=fontsize,color='black')
    
    if mirror_ang > 90:
        txt='  mirror angle = %.2f$^\circ$' % (mirror_ang-90)
    else:
        txt='  mirror angle = %.2f$^\circ$' % mirror_ang        
    ax.text(x[0],y[0],txt,ha='left',va='top',fontsize=fontsize,color='black')

    # dashed line to show mirror position
    x=[mirror_pos[0],mirror_pos[0]]
    y=[mirror_pos[1],0]
    ax.plot(x,y,color='black',linestyle='dashed')

    x=mirror_pos[0]
    y=mirror_pos[1]
    txt=' flat mirror at\n%.2fm, %.2fm' % (abs(mirror_pos[0]),mirror_pos[1])
    ax.text(x,y,txt,color='black',ha='right',va='bottom',fontsize=fontsize)

    # total optical path
    Dcal_mirror=sqrt( (cal_pos[0]-mirror_pos[0])**2 + (cal_pos[1]-mirror_pos[1])**2  )
    Dqubic_mirror=sqrt( (mirror_pos[0]-window_pos[0])**2 + (mirror_pos[1]-window_pos[1])**2 )
    optical_dist=Dcal_mirror+Dqubic_mirror

    # edge line
    x=[0,0]
    y=[0,Hsource]
    #ax.plot(x,y,color='blue')

    # bottom line
    x=[1-xlength,1]
    y=[0,0]
    ax.plot(x,y,color='black',linewidth=0.5,linestyle='dotted')

    # horizontal line at window height to compare to the drawing by Didier
    x=[1-xlength,1]
    y=[window_pos[1],window_pos[1]]
    ax.plot(x,y,color='black',linewidth=0.5,linestyle='dotted')

    # draw a box outline of the QUBIC cryostat
    top_centre=[window_pos[0]-win_offset*sin(pointing_ang_rad),window_pos[1]-win_offset*cos(pointing_ang_rad)]
    bot_centre=[top_centre[0]+cryo_height*cos(pointing_ang_rad),top_centre[1]-cryo_height*sin(pointing_ang_rad)]
    x=[bot_centre[0],top_centre[0]]
    y=[bot_centre[1],top_centre[1]]
    ax.plot(x,y,color='black',linewidth=1,linestyle='dotted')
    #ax.plot(top_centre[0],top_centre[1],marker='X',color='black')

    cryo_corner1=(top_centre[0]-0.5*cryo_width*sin(pointing_ang_rad),top_centre[1]-0.5*cryo_width*cos(pointing_ang_rad))
    cryo_corner2=(top_centre[0]+0.5*cryo_width*sin(pointing_ang_rad),top_centre[1]+0.5*cryo_width*cos(pointing_ang_rad))
    cryo_corner3=(bot_centre[0]+0.5*cryo_width*sin(pointing_ang_rad),bot_centre[1]+0.5*cryo_width*cos(pointing_ang_rad))
    cryo_corner4=(bot_centre[0]-0.5*cryo_width*sin(pointing_ang_rad),bot_centre[1]-0.5*cryo_width*cos(pointing_ang_rad))
    pts=[cryo_corner1,cryo_corner2,cryo_corner3,cryo_corner4,cryo_corner1]
    x=[val[0] for val in pts]
    y=[val[1] for val in pts]
    ax.plot(x,y,color='black',linestyle='dotted')
    
    # plot extents
    ax.set_xlim(1-xlength,1)
    ax.set_ylim(-1,-1+ylength)
    # equal aspect ratio
    #ax.set_aspect('equal', 'datalim')
    ax.set_aspect('equal')

    # labels
    txt='Total optical path length = %.2fm' % optical_dist
    ax.text(0.5,0.1,txt,va='top',ha='center',fontsize=fontsize,transform=ax.transAxes)

    ax.text(0.5,1.01,'QUBIC Calibration Source layout',ha='center',va='bottom',fontsize=fontsize,transform=ax.transAxes)
    ax.tick_params(labelsize=fontsize)
    
    pdfname='/home/work/qubic/pix/calsource_setup/QUBIC_calibration_source_layout_sketch.pdf'
    fig.savefig(pdfname,format='pdf',dpi=300,bbox_inches='tight',transparent=transparency)
    return
