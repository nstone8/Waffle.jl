module Waffle
using Tessen, Delica, Unitful, Ahmes, Statistics, Primes
#go get some internal Tessen stuff
import Tessen:HatchLine, pointalong, intersections

export createconfig, scaffold, psweep

"""
```julia
createconfig([filename])
```
Write an example config file to `filename`. If `filename` is omitted, the file
will be written to `"config.jl"`.

# Configuration Parameters
- lscaf: scaffold length
- wscaf: scaffold width
- dfield: calibrated FOV for the objective being used
- hbottom: distance from the bottom of the post to the bottom of the beams
- htop: distance from the bottom of the beams to the top of the posts
- chamferbottom: angle at which the bottoms of the posts and beams should be chamfered
- chamfertop: angle at which the topss of the posts and beams should be chamfered
- cutangle: angles at which blocks should be cut to avoid shadowing
- overlap: amount that neighboring blocks should overlap to ensure they are connected
- dslice: slicing distance
- wbumper: width of the bumpers
- fillet: fillet radius on xz and yz crossections of the posts and bumpers
- wpost: post width
- hbeam: beam height
- wbeam: beam width
- lbeammax: maximum total beam length
- maxseglength: maximum length of a single beam segment
- keygap: closest point between the two halves of a beam before they are closed with a keystone
- dhammockhatch: hammock hatching distance
- dhatch: hatching distance for posts, beams and bumpers
- laserpower: laser power for posts, beams and bumpers
- scanspeed: scan speed for posts, beams and bumpers
- stagespeed: max stage speed
- interfacepos: how far into the substrate writing should begin
- hamscanspeed: scan speed for hammocks
- hamlaserpower: laser power for hammocks
"""
function createconfig(filename="config.jl")
    config = """Dict(
        :lscaf => 6u"mm",
        :wscaf => 4u"mm",
        :dfield => 1750u"µm",
        :hbottom => 30u"µm",
        :htop => 50u"µm",
        :chamferbottom => pi/6,
        :chamfertop => pi/6,
        :cutangle => pi/6,
        :overlap => 10u"µm",
        :dslice => 1u"µm",
        :wbumper => 200u"µm",
        :fillet => 30u"µm",
        :wpost => 100u"µm",
        :hbeam => 10u"µm",
        :wbeam => 25u"µm",
        :lbeammax => 120u"µm",
        :maxseglength => 30u"µm",
        :keygap => 20u"µm",
        :dhammockhatch => 1u"µm",
        :dhatch => 300u"nm",
        :laserpower => 100u"mW",
        :scanspeed => 100000u"µm/s",
        :stagespeed => 100u"µm/s",
        :interfacepos => 10u"µm",
        :hamscanspeed => 10000u"µm/s",
        :hamlaserpower => 100u"mW"
    )
    """
    open(filename,"w") do io
        print(io,config)
    end
end

function bumperedgecoords(;kwargs...)
    #copy out the variables we want
    hbottom = kwargs[:hbottom]
    htop = kwargs[:htop]
    chamfertop = kwargs[:chamfertop]
    chamferbottom = kwargs[:chamferbottom]
    wbumper = kwargs[:wbumper]
    fillet = kwargs[:fillet]
    #figuring out how to fillet everything in crossection seems hard enough that I'm just going
    #to use Tessen hatching machinery to do it for me.
    #this function will return a closure which gives the bottom and top coordinates of a bumper
    #as a function of z position

    #draw our crossection on the xz plane
    verts=[
        [-wbumper/2,hbottom],
        [-wbumper/2 + htop*tan(chamfertop),hbottom+htop],
        [wbumper/2 - htop*tan(chamfertop),hbottom+htop],
        [wbumper/2,hbottom],
        [wbumper/2 - hbottom*tan(chamferbottom),0u"µm"],
        [-wbumper/2 + hbottom*tan(chamferbottom),0u"µm"]
    ]
    nofillet = polycontour(verts)
    filleted = polycontour(verts,fillet)
    function(zcoord)
        #place a hatchline off of the contour to the left and at our desired height
        hl=HatchLine([-wbumper,zcoord],[1,0]u"µm")
        #our first coordinate is the first intersection of hl with nofillet
        firstpara = sort(intersections(nofillet,hl))[1]
        #our second coordinate is the second intersection of hl with filleted
        secondpara = sort(intersections(filleted,hl))[2]
        map([firstpara,secondpara]) do p
            #convert to real world coordinates in µm
            coords = pointalong(hl,p)
            #we're interested in the x coordinate, add the dimensions on
            coords[1] * u"µm"
        end
    end
end

#helper function to build a series of segments to get us from startx to endx
#we will assume startx and endx do not include overlap on either side
function bumpersegtrain(startx,endx,startz,endz,bec;kwargs...)
    #to connect the left and right sides, we will use a bunch of identical segments cut to
    #have overhangs on the left side and placed to overlap with each other
    #get the amount of length that should be devoted to overhangs
    lseg = sqrt(kwargs[:dfield]^2 - kwargs[:wbumper]^2)
    olength = 2*(endz-startz)*tan(kwargs[:cutangle])
    #the length of every slice should be constant, the maximum progress we can make per slice
    #is set by
    lslicemax = lseg - (olength/2) - kwargs[:overlap]
    #the distance we have to cover (including overlapping at either end) is...
    reqdist = endx -  startx + 2*sign(endx-startx)*kwargs[:overlap]

    numsegs = ceil(Int,abs(reqdist/lslicemax))
    distperseg = reqdist/numsegs

    #build a block representing one of these segments with the middle of the bottom of the `source`
    #edge at the local origin

    segslices = map(range(start=startz,step=kwargs[:dslice],stop=endz)) do z
        ycoords = bec(z)
        #we should slant the opposite direction when moving backwards
        sliceoffset = -sign(reqdist)*z*tan(kwargs[:cutangle])
        #draw the edges
        e = [LineEdge([sliceoffset,ycoords[1]],
                      [sliceoffset+distperseg+sign(distperseg)*kwargs[:overlap],ycoords[1]]),
             LineEdge([sliceoffset+distperseg+sign(distperseg)*kwargs[:overlap],ycoords[1]],
                      [sliceoffset+distperseg+sign(distperseg)*kwargs[:overlap],ycoords[2]]),
             LineEdge([sliceoffset+distperseg+sign(distperseg)*kwargs[:overlap],ycoords[2]],
                      [sliceoffset,ycoords[2]]),
             LineEdge([sliceoffset,ycoords[2]],
                      [sliceoffset,ycoords[1]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    #place the first such block overlapping with lbblock
    firstseg = Block(segslices...,origin=[startx - sign(distperseg)*kwargs[:overlap],
                                          0u"µm",startz])
    #we actually want the objective to be centered on the segment during the print, so we need to
    #move the local origin. we will do this by sliding the geometry to the left while using
    #preserveframe, and then slide the resulting object back
    #the -startz translation in the preserveframe step moves the geometry to z=0 in the local
    #frame, the origin of which is at startz in the enclosing frame
    fbsdisp = sign(distperseg)*olength+(distperseg/2)
    blocks = [translate(translate(firstseg,[-fbsdisp,0u"µm",-startz],preserveframe=true),
                        [fbsdisp,0u"µm"])]
    #build the rest by copying
    for _ in 2:numsegs
        push!(blocks,translate(blocks[end],[distperseg,0u"µm"]))
    end
    return blocks
end

function bumper(;kwargs...)
    #The bumper is going to be way longer than `dfield`. The longest segment with width `w`
    #we can write at a time is
    lseg = sqrt(kwargs[:dfield]^2 - kwargs[:wbumper]^2)
    #this includes all of the space we have for overlaps

    #get our width profile as a function of z coordinate
    bec=bumperedgecoords(;kwargs...)
    #we first need to build the bottom of the leftmost end of the bumper. the width of this block will come
    #from `bec`, the left face will be completely rounded off.
    #we would like the right face of the block to taper down to the substrate at `cutangle`
    #we would like this taper to hit the substrate at x=lseg/2, and for the rounded arc of the
    #widest slice (which will have a radius of `wbumper` to be tangent to x=-lseg/2

    #we will build the bumper up to a z height of (hbottom+htop+overlap)/2 on the first pass, and
    #then glue on the top coming back
    zmid=(kwargs[:hbottom]+kwargs[:htop]+kwargs[:overlap])/2
    lbslices = map(range(start=0u"µm",step=kwargs[:dslice],stop=zmid)) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom right corner, this should be lseg/2
        #for the bottom slice
        xbr = (lseg/2) - z*tan(kwargs[:cutangle])
        #get the center of the arc making up the left face
        cleftarc = [-lseg/2 + kwargs[:wbumper]/2, mean(ycoords)]
        #draw the edges
        e = [ArcEdge(cleftarc,wslice/2,pi/2,3pi/2),
             LineEdge([cleftarc[1],ycoords[1]],[xbr,ycoords[1]]),
             LineEdge([xbr,ycoords[1]],[xbr,ycoords[2]]),
             LineEdge([xbr,ycoords[2]],[cleftarc[1],ycoords[2]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    #we want the origin of this block to be at [(-lscaf+lseg)/2,0] so that the whole object is
    #centered on the origin
    lbblock=Block(lbslices...,origin=[(-kwargs[:lscaf]+lseg)/2,0u"µm",0u"µm"])

    #we will now do the bottom of the rightmost block. This will be the mirror image of lbblock
    #but with the cut going the opposite way (i.e. the topmost slice is longest)
    rbslices = map(range(start=0u"µm",step=kwargs[:dslice],stop=zmid)) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom left corner, this should be -lseg/2
        #for the topmost slice
        xbl = (-lseg/2) + (zmid-z)*tan(kwargs[:cutangle])
        crightarc = [lseg/2 - kwargs[:wbumper]/2, mean(ycoords)]
        #draw the edges
        e = [ArcEdge(crightarc,wslice/2,3pi/2,pi/2),
             LineEdge([crightarc[1],ycoords[2]],[xbl,ycoords[2]]),
             LineEdge([xbl,ycoords[2]],[xbl,ycoords[1]]),
             LineEdge([xbl,ycoords[1]],[crightarc[1],ycoords[1]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    #we want the origin of this block at [(lscaf-lseg)/2,0]
    rbblock=Block(rbslices...,origin=[(kwargs[:lscaf]-lseg)/2,0u"µm",0u"µm"])

    #to connect the left and right sides, we will use a bunch of identical segments cut to
    #have overhangs on the left side and placed to overlap with each other
    midbottomblocks = bumpersegtrain(-kwargs[:lscaf]/2 + lseg,
                                     kwargs[:lscaf]/2 - lseg + zmid*tan(kwargs[:cutangle]),
                                     0u"µm",zmid,bec;kwargs...)
    #now we have to build the top. We'll do similar to what we did for the bottom, but we'll
    #make the end pieces sorter so the seams between blocks on the top and bottom are less
    #likely to line up

    zmid2=(kwargs[:hbottom]+kwargs[:htop]-kwargs[:overlap])/2
    rtslices = map(range(start=zmid2,step=kwargs[:dslice],stop=kwargs[:hbottom]+kwargs[:htop])) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom left corner, this should be -lseg/4
        #for the bottom slice
        xbl = (-lseg/4) + (z)*tan(kwargs[:cutangle])
        crightarc2 = [lseg/4 - kwargs[:wbumper]/2, mean(ycoords)]
        #draw the edges
        e = [ArcEdge(crightarc2,wslice/2,3pi/2,pi/2),
             LineEdge([crightarc2[1],ycoords[2]],[xbl,ycoords[2]]),
             LineEdge([xbl,ycoords[2]],[xbl,ycoords[1]]),
             LineEdge([xbl,ycoords[1]],[crightarc2[1],ycoords[1]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    rtblock=Block(rtslices...,origin=[kwargs[:lscaf]/2-lseg/4,0u"µm",zmid2])
    #given that we want the origin of the coordinate system at zmid2, we need to translate
    #the geometry in rtblock to start at z=0 in the local frame
    rtblock=translate(rtblock,[0u"µm",0u"µm",-zmid2],preserveframe=true)

    ltslices = map(range(start=zmid2,step=kwargs[:dslice],stop=kwargs[:hbottom]+kwargs[:htop])) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom right corner, this should be lseg/4
        #for the topmost slice
        xbr = (lseg/4) + (z-kwargs[:htop])*tan(kwargs[:cutangle])
        cleftarc2 = [-lseg/4 + kwargs[:wbumper]/2, mean(ycoords)]
        #draw the edges
        e = [ArcEdge(cleftarc2,wslice/2,pi/2,3pi/2),
             LineEdge([cleftarc2[1],ycoords[1]],[xbr,ycoords[1]]),
             LineEdge([xbr,ycoords[1]],[xbr,ycoords[2]]),
             LineEdge([xbr,ycoords[2]],[cleftarc2[1],ycoords[2]])]
        slice = Slice([Contour(e)])
        z => slice
    end
    ltblock=Block(ltslices...,origin=[-kwargs[:lscaf]/2+lseg/4,0u"µm",zmid2])
    #given that we want the origin of the coordinate system at zmid2, we need to translate
    #the geometry in ltblock to start at z=0 in the local frame
    ltblock=translate(ltblock,[0u"µm",0u"µm",-zmid2],preserveframe=true)

    #build the middle segments in the top
    midtopblocks = bumpersegtrain(kwargs[:lscaf]/2 - lseg/2,
                                  -kwargs[:lscaf]/2 + lseg/2 -
                                      (kwargs[:htop]-zmid2)*tan(kwargs[:cutangle]),
                                  zmid2,kwargs[:hbottom]+kwargs[:htop],bec;kwargs...)
    return SuperBlock(lbblock,midbottomblocks...,rbblock,rtblock,midtopblocks...,ltblock)
end

#get one slice of a post
function postslice(scale,pctfillet;kwargs...)
    #define some vertices that we can use to get the rest via rotation
    firstverts = scale .* [[-kwargs[:wbeam]/2,kwargs[:wpost]/2],
                           [kwargs[:wbeam]/2,kwargs[:wpost]/2]]
    restverts = map((-pi/2):(-pi/2):(-3pi/2)) do r
        map(firstverts) do v
            Tessen.zrotate(v,r)
        end
    end
    verts=vcat(firstverts,restverts...)
    #calculate the maximum possible fillet
    vlongside = verts[3] - verts[2]
    phi = atan(/(reverse(vlongside)...)) |> abs
    theta = pi - phi
    rmax = scale*kwargs[:wbeam] * tan(theta/2)/2
    Slice([polycontour(verts,rmax*pctfillet)])
end

function post(;kwargs...)
    #get all the z slices that we will need
    unassignedz = collect(range(0u"µm",kwargs[:hbottom]+kwargs[:htop],step=kwargs[:dslice]))
    bottomz = filter(unassignedz) do uz
        uz < kwargs[:hbottom]
    end
    filter!(unassignedz) do uz
        !(uz in bottomz)
    end
    midz = filter(unassignedz) do uz
        uz <= kwargs[:hbottom] + kwargs[:hbeam]
    end

    filter!(unassignedz) do uz
        !(uz in midz)
    end
    
    topz = unassignedz
    #bottom will go from fully rounded at the bottom to no rounding at z=hbottom
    bottomslices = map(bottomz) do z
        #the amount of undercut
        ucut = (kwargs[:hbottom]-z)*tan(kwargs[:chamferbottom])
        #in terms of percent of wpost
        scale = (kwargs[:wpost]-ucut)/kwargs[:wpost]
        pctround = (kwargs[:hbottom] - z)/kwargs[:hbottom]
        z => postslice(scale,pctround;kwargs...)
    end
    #midslices has no rounding so the beams will snap on easily
    midslices = map(midz) do z
        #the amount of undercut
        ucut = (z-kwargs[:hbottom])*tan(kwargs[:chamfertop])
        #in terms of percent of wpost
        scale = (kwargs[:wpost]-ucut)/kwargs[:wpost]
        z => postslice(scale,0;kwargs...)
    end
    #topslices round smoothly up to 1
    topslices = map(topz) do z
        #the amount of undercut
        ucut = (z-kwargs[:hbottom])*tan(kwargs[:chamfertop])
        #in terms of percent of wpost
        scale = (kwargs[:wpost]-ucut)/kwargs[:wpost]
        pctround = (z - kwargs[:hbottom] - kwargs[:hbeam]) / (kwargs[:htop] - kwargs[:hbeam])
        z => postslice(scale,pctround;kwargs...)
    end
    Block(bottomslices...,midslices...,topslices...)
end

#build a segmented beam with nsegs segments centered on y=0.
#startx and stopx should be the position of the edge we are bonding to at the z coordinate
#corresponding to the bottom of the beam (i.e. this function will build in overlap)
#gonna go overboard and make this a struct so we can define rotation and translation
struct Beam
    leftsegs
    rightsegs
    keystone
end

function Beam(nsegs::Int,startx,stopx;kwargs...)
    #nsegs needs to be odd
    @assert isodd(nsegs)
    #we want the keystone to be as short as possible. Point of closest approach between the
    #two halves is set by `keygap`
    #length of the keystone measured at the beam midline
    lkey = kwargs[:keygap] + kwargs[:hbeam]*tan(kwargs[:cutangle])
    #the first segment needs to be chamfered according to `chamfertop`, we will do the rest at
    #cutangle
    distperseg = (stopx - startx - lkey)/(nsegs-1)
    lseg = distperseg + kwargs[:overlap]
    leftsegs = map(1:((nsegs-1)/2)) do i
        #the first segment needs to be a little longer and chamfered differently
        segpos = startx + distperseg*((2i-1)/2)
        seg = if i==1
            box(lseg + (kwargs[:overlap]/2),kwargs[:wbeam],kwargs[:hbeam],kwargs[:dslice],
                chamfer =[-kwargs[:chamfertop] kwargs[:cutangle]
                          kwargs[:cutangle]    kwargs[:cutangle]])            
        else
            box(lseg,kwargs[:wbeam],kwargs[:hbeam],kwargs[:dslice],
                chamfer =[-kwargs[:cutangle]   kwargs[:cutangle]
                          kwargs[:cutangle]    kwargs[:cutangle]])            
        end
        #seg is currently centered at [0,0,0]. Move it into position (use preserveframe so we don't
        #move the stage
        seg=translate(seg,[segpos,0u"µm",kwargs[:hbottom]+kwargs[:hbeam]/2],preserveframe=true)
        if i==1
            #scrootch a little bit back to make the overlap right
            seg=translate(seg,[-kwargs[:overlap]/4,0u"µm",0u"µm"],preserveframe=true)
        end
        return seg
    end
    #the center of the beam in xy plane
    cbeam = [(startx+stopx)/2,0u"µm"]
    #we can make the right side of the bridge by rotating leftsegs around cbeam
    rightsegs = map(leftsegs) do ls
        rotate(ls,pi,cbeam,preserveframe=true)
    end
    #the keystone is cut differently
    keystone = box(lkey,kwargs[:wbeam],kwargs[:hbeam],kwargs[:dslice],
                   chamfer =[-kwargs[:cutangle]   -kwargs[:cutangle]
                             kwargs[:cutangle]    kwargs[:cutangle]])
    #move it into position
    keystone = translate(keystone,
                         vcat(cbeam,kwargs[:hbottom]+kwargs[:hbeam]/2),
                         preserveframe=true)
    #return a namedtuple so we can be fancy about how we write these
    Beam(leftsegs,rightsegs,keystone)
end

function Tessen.translate(b::Beam,args...;kwargs...)
    Beam([translate(x,args...;kwargs...) for x in b.leftsegs],
         [translate(x,args...;kwargs...) for x in b.rightsegs],
         translate(b.keystone,args...;kwargs...))
end

function Tessen.rotate(b::Beam,args...;kwargs...)
    Beam([rotate(x,args...;kwargs...) for x in b.leftsegs],
         [rotate(x,args...;kwargs...) for x in b.rightsegs],
         rotate(b.keystone,args...;kwargs...))
end

#draw a hammock, it will be centered on (0,0) if topflat and bottomflat are both false
#those arguments will add corners so it mates flush with a bumper
function hammock(lbeamx,lbeamy,topflat,bottomflat;kwargs...)
    #the coordinates of the vertices which make up the top right cut off corner are
    firstverts = [[lbeamx/2,(lbeamy + kwargs[:wpost] - kwargs[:wbeam])/2],
                  [(lbeamx + kwargs[:wpost] - kwargs[:wbeam])/2,lbeamy/2]]
    #get all the vertices by mirroring
    cornerfactors=[[1,1],[1,-1],[-1,-1],[-1,1]]
    invertverts=[false,true,false,true]
    allverts = map(zip(cornerfactors,invertverts)) do (cf,iv)
        newverts = [cf.*v for v in firstverts]
        iv ? reverse(newverts) : newverts
    end
    itoflat=[]
    if topflat
        #flatten out the top two corners
        push!(itoflat,1,4)
    end
    if bottomflat
        push!(itoflat,2,3)
    end
    for i in itoflat
        oldverts = allverts[i]
        absmat = abs.(hcat(oldverts...))
        #take the maximum coordinate in x, and the minimum in y
        xmag=maximum(absmat[1,:])
        ymag=minimum(absmat[2,:])
        #get the right sign
        allverts[i] = [[sign(oldverts[1][1])*xmag,sign(oldverts[1][2])*ymag]]
    end
    #flatten to a vector of vectors
    flatverts = vcat(allverts...)
    #scale coordinates to give the correct overlap on the beams
    xscale = (kwargs[:overlap]+(lbeamx + kwargs[:wpost] - kwargs[:wbeam])/2)/
        ((lbeamx + kwargs[:wpost] - kwargs[:wbeam])/2)
    yscale = (kwargs[:overlap]+(lbeamy + kwargs[:wpost] - kwargs[:wbeam])/2)/
        ((lbeamy + kwargs[:wpost] - kwargs[:wbeam])/2)
    scaledverts = map(flatverts) do fv
        fv .* [xscale,yscale]
    end
    outline = Slice([polycontour(scaledverts)])
    hatchedslices = map([pi/4,3pi/4]) do hd
        (kwargs[:hbottom]+kwargs[:hbeam]/2) =>
            hatch(outline,dhatch=kwargs[:dhammockhatch],hatchdir=hd)
    end
    Block(hatchedslices...)
end

#build a kernel with a nx x ny array of posts with pitch of px and py in each dimension
#lbeams and bbeams controls wheter beams should be drawn connecting to geometry to the left
#or below the current kernel
function kernel(nx,ny,px,py,nseg;lbeams,bbeams,toprow,kwargs...)
    #calculate beam lengths
    lbeamx = px - kwargs[:wpost]
    lbeamy = py - kwargs[:wpost]
    #with the optional beams to the left and bottom, the total size of the array is...
    kernelsize = [nx*(kwargs[:wpost] + lbeamx) + kwargs[:overlap], #one dangling beam
                  ny*(kwargs[:wpost] + lbeamy) + lbeamy + 2*kwargs[:overlap]] #two danglers
    #if we think of the upper left corner of this maximal array as (0,0) the post coordinates are...
    #increasing y is up/north
    topleftpostcoords = [[i*(lbeamx + kwargs[:wpost]) - kwargs[:wpost]/2 + kwargs[:overlap],
                          -j*(lbeamy + kwargs[:wpost]) +
                              kwargs[:wpost]/2 -kwargs[:overlap]] for i in 1:nx, j in 1:ny]
    #however, we want the objective at the middle of this maximal array, so let's move all
    #these points to be centered around the center of the 'maximal' array
    postcoords = map(topleftpostcoords) do tlpc
        tlpc - [1,-1].*kernelsize/2
    end
    #little helper function
    function onepost(postcenter,left,bottom,toprow)
        thispost = translate(post(;kwargs...),postcenter,preserveframe=true)
        #beam(nsegs::Int,startx,stopx;kwargs...)
        #first build 'raw' beams around (0,0) and then translate them
        rawleftbeam = Beam(nseg,-kwargs[:wpost]/2 - lbeamx,
                           -kwargs[:wpost]/2;kwargs...)
        #topbeam will also require rotation (beam() makes horizontal beams)
        rawtopbeam = Beam(nseg,kwargs[:wpost]/2,kwargs[:wpost]/2 + lbeamy;kwargs...)
        leftbeam = translate(rawleftbeam,postcenter,preserveframe=true)
        topbeam = translate(rotate(rawtopbeam,pi/2,preserveframe=true),
                            postcenter,preserveframe=true)
        allbeams=[topbeam]
        hammocks=Block{HatchedSlice}[]
        if left
            push!(allbeams,leftbeam)
            ht = hammock(lbeamx,lbeamy,toprow,false;kwargs...)
            push!(hammocks,translate(ht,postcenter + [-(kwargs[:wpost] + lbeamx)/2,
                                                      (kwargs[:wpost] + lbeamy)/2],
                                     preserveframe=true))
            #also need one underneath us if bottom
            if bottom
                #i.e. if left and bottom, this implies we are on the bottom row
                hb=hammock(lbeamx,lbeamy,false,true;kwargs...)
                push!(hammocks,translate(hb,postcenter + [-(kwargs[:wpost] + lbeamx)/2,
                                                          -(kwargs[:wpost] + lbeamy)/2],
                                         preserveframe=true))
            end
        end
        if bottom
            #make the bottom beam by rotating topbeam
            push!(allbeams,rotate(topbeam,pi,postcenter,preserveframe=true))
        end
        #make a namedtuple so we can arrange things nicely in the enclosing scope
        (post=thispost,beams=allbeams,hammocks=hammocks)
    end
    needleft = [(i==1) ? lbeams : true for i in 1:nx, j in 1:ny]
    needbottom = [(j==ny) ? bbeams : false for i in 1:nx, j in 1:ny]
    intoprow = [(j==1) ? toprow : false for i in 1:nx, j in 1:ny]
    units = map(zip(postcoords,needleft,needbottom,intoprow)) do (pc,nl,nb,tr)
        onepost(pc,nl,nb,tr)
    end
    uvec = reshape(units,:)
    #get a vector where each entry is a 'half' beam
    beamvec = vcat(map(uvec) do uv
                       vcat(map(uv.beams) do b
                                [b.leftsegs, b.rightsegs]
                            end...)
                   end...)
    keystones = collect(b.keystone for u in uvec for b in u.beams)
    posts = [u.post for u in uvec]
    hammocks = vcat((u.hammocks for u in uvec)...)
    #now we arrange how we want things to be printed
    #all the posts should be written in parallel
    postblock=merge(posts...)
    #now we want each 'matching' segment of the beams to be written in parallel
    #i.e. the first leftseg and the first rightseg for each beam should be written together,
    #then the second, etc. all of the keystones last
    segblocks = map(zip(beamvec...)) do thesesegs
        merge(thesesegs...)
    end
    keyblock = merge(keystones...)
    #so now we print the posts first, then the segs, then the keystones, and then the hammocks
    (posts=postblock,support=SuperBlock(segblocks...,keyblock), hammocks=SuperBlock(hammocks...))
end

"""
```julia
scaffold(scaffolddir[, configfilename])
scaffold(scaffolddir,configdict)
```
Build a directory of .gwl files to build a scaffold in a directory at `scaffolddir`. The
scaffold's geometrical parameters can be provided as a `Dict` or read from `configfilename`.
If `configfilename` is omitted, the parameters will be read from `"config.jl"`
"""
function scaffold end

function scaffold(scaffolddir,kwargs::Dict)
    #make a folder to hold the files for this scaffold
    mkdir(scaffolddir)
    #move down in there do compile the geometry
    compgeom = cd(scaffolddir) do
        #make a folder for our pre-compiled scaffold parts
        mkdir("scripts")
        #assemble and pre-compile the bumpers
        b=bumper(;kwargs...)
        #this bumper is currently at the origin, translate it up into position
        btop=translate(b,[0u"µm",(kwargs[:wscaf]-kwargs[:wbumper])/2,0u"µm"])
        #get the bottom by rotation
        bbottom = rotate(btop,pi)
        @info "compiling bumpers"
        tophatched=hatch(btop,dhatch=kwargs[:dhatch],bottomdir=pi/4)
        bottomhatched=hatch(bbottom,dhatch=kwargs[:dhatch],bottomdir=pi/4)
        #printing the 'bottom' bumper first reduces stage movement
        bumpers = CompiledGeometry(joinpath("scripts","bumpers.gwl"),bottomhatched,tophatched;
                                   laserpower=kwargs[:laserpower],
                                   scanspeed=kwargs[:scanspeed])
        #determine our kernel parameters
        #maximum possible pitch
        pmax = kwargs[:wpost] + kwargs[:lbeammax]
        #calculate our grid dimensions
        
        ny = ceil(Int,
                  #distance from the top of the top row of posts to the bottom bumper
                  (kwargs[:wscaf] - 2*kwargs[:wbumper] - kwargs[:lbeammax])/
                      pmax)
        
        nx = #start with one post
            1 + 
            ceil(Int,
                 #distance from the end of the leftmost post to the end of the scaffold
                 #we have to subtract off wbumper to account for the rounded ends of the
                 #bumper
                 (kwargs[:lscaf] - kwargs[:wpost] - kwargs[:wbumper])/
                     pmax)

        #we don't want nx or ny to be prime
        (nx,ny) = map([nx,ny]) do ni
            isprime(ni) ? ni+1 : ni
        end
        
        #get the actual beam lengths
        nbeamx = nx - 1
        nbeamy = ny + 1
        lbeamx = (kwargs[:lscaf] - nx*kwargs[:wpost] - kwargs[:wbumper])/nbeamx
        lbeamy = (kwargs[:wscaf]-2*kwargs[:wbumper] - ny*kwargs[:wpost])/nbeamy
        @assert lbeamx <= kwargs[:lbeammax]
        @assert lbeamy <= kwargs[:lbeammax]
        #calculate actual pitch
        px = lbeamx + kwargs[:wpost]
        py = lbeamy + kwargs[:wpost]
        #figure out the biggest kernel which fits neatly into our grid
        #define a helper function so we can get all multiples of our nx and ny
        function allmultiples(x)
	    filter(1:x) do y
	        (x%y) == 0
	    end
        end

        #what is the biggest 'kernel' of posts that 1) fits evenly into the desired space
        #and can fit inside a circle with diameter `dfield`
        #we will assume each 'kernel' includes connecting bridges to the left, top and bottom.
        #the left bridges will not be written on the left edge and the bottom bridges will only
        #actually be written on the bottom row

        allnumcombos = [(nxi,nyi) for nxi in allmultiples(nx), nyi in allmultiples(ny)]
        possiblenumcombos = filter(allnumcombos) do (nxi,nyi)
	    sizex = nxi*px + kwargs[:overlap] #one 'dangling' beam
	    sizey = nyi*py + lbeamy + 2*kwargs[:overlap] #two 'dangling' beams
	    sqrt(sizex^2 + sizey^2) < kwargs[:dfield] #does this fit
        end
        #get the option which contains the most posts
        (_,maxcomboindex) = map(possiblenumcombos) do (nx,ny)
	    nx*ny
        end |> findmax
        #assign this max size to variables for convenience
        (knx,kny) = possiblenumcombos[maxcomboindex]
        #pre-compile all the different kernel variations assuming we have more than one FOV in every
        #direction
        #first we need to know the number of segments we're going to want in our beams
        longestbeam = max(lbeamx,lbeamy)
        #this is the definition of lkey copy-pasted from beam()
        lkey = kwargs[:keygap] + kwargs[:hbeam]*tan(kwargs[:cutangle])
        halfseglength = (longestbeam - lkey)/2
        #did this using half the beam to guarentee this would be a even number
        nseg = 2*(ceil(Int,halfseglength/kwargs[:maxseglength]))
        #we count the keystone as a segment above
        nseg += 1
        #going from left to right, we care if a kernel is on the leftmost edge or not as kernels
        #on the leftmost edge don't need left beams
        xvariants = Dict(:leftedge => (lbeams = false,),:notleftedge => (lbeams = true,))
        #from top to bottom, we care if we're on the top, in the middle, or on the bottom
        #we will also build the case where we're on the top and bottom (i.e. the scaffold is one
        #FOV in width for completeness
        yvariants = Dict(
            :top => (toprow=true,bbeams=false),
            :middle => (toprow=false,bbeams=false),
            :bottom => (toprow=false,bbeams=true),
            :topbottom => (toprow=true,bbeams=true)
        )
        #build all the combinations
        kvariants = Dict()
        #all posts are the same
        postblock=kernel(knx,kny,px,py,nseg;
                         xvariants[:notleftedge]...,yvariants[:middle]...,kwargs...).posts
        @info "compiling posts"
        hatchedposts=hatch(postblock,dhatch=kwargs[:dhatch],bottomdir=pi/4)
        posts=CompiledGeometry(joinpath("scripts","posts.gwl"),SuperBlock(hatchedposts);
                               laserpower=kwargs[:laserpower],
                               scanspeed=kwargs[:scanspeed])
        
        for xv in keys(xvariants)
            kvariants[xv] = Dict()
            for yv in keys(yvariants)
                thiskernel = kernel(knx,kny,px,py,nseg;
                                    xvariants[xv]...,yvariants[yv]...,kwargs...)
                #hatch it
                @info "compiling $(join([xv,yv])) kernel"
                thissupport = hatch(thiskernel.support;dhatch=kwargs[:dhatch],bottomdir=pi/4)
                thishammocks = hatch(thiskernel.hammocks;dhatch=kwargs[:dhatch],bottomdir=pi/4)
                tksupport = CompiledGeometry(joinpath("scripts",join([xv,yv,"_structure.gwl"])),
                                             thissupport;laserpower=kwargs[:laserpower],
                                             scanspeed=kwargs[:scanspeed])
                
                tkhammocks = CompiledGeometry(joinpath("scripts",join([xv,yv,"_hammocks.gwl"])),
                                              thishammocks;laserpower=kwargs[:hamlaserpower],
                                              scanspeed=kwargs[:hamscanspeed])
                
                kvariants[xv][yv]=Dict(:support => tksupport,:hammocks => tkhammocks)
            end
        end
        #place all the kernels and hammocks
        #the kernels were built so that they are centered on a 'maximal' kernel
        #the coordinates of the top left corner of the top left post in the top left kernel relative
        #to the top left corner of the scaffold is...
        cornertofirstpost = [kwargs[:wbumper]/2,-kwargs[:wbumper] - lbeamy]
        #the distance between the top left corner of the first post and the origin of its kernel is
        firstposttokc = [(knx*px)/2-lbeamx,(-kny*py-lbeamy)/2 + lbeamy]
        #therefore, the top left corner of the scaffold to the origin of the first kernel is
        cornertokc = cornertofirstpost + firstposttokc
        numkernelx = nx/knx
        numkernely = ny/kny
        #the coords of all of the kernel centers relative to the top left corner
        cornertoallcenters = [cornertokc + [knx*px*i,-kny*py*j] for
                                  i in 0:(numkernelx-1), j in 0:(numkernely-1)]
        #coords of all of the kernel origins relative to the center of the scaffold
        relcenters = map(cornertoallcenters) do ctac;
            ctac - [kwargs[:lscaf]/2, -kwargs[:wscaf]/2]
        end
        #place our kernel objects in the right positions
        kernelindices = [[i,j] for i in 1:numkernelx, j in 1:numkernely]
        kernelmat = map(zip(kernelindices,relcenters)) do ((i,j),kcenter)
            #add the z coordinate to kcenter
            push!(kcenter,0u"µm")
            #figure out what kernel variant we want
            xv = (i==1) ? :leftedge : :notleftedge
            yv = if numkernely==1
                :topbottom
            elseif j==1
                :top
            elseif j==numkernely
                :bottom
            else
                :middle
            end
            tk=kvariants[xv][yv]
            #need to translate both of them to kcenter
            tposts = translate(posts,kcenter)
            tsupport=translate(tk[:support],kcenter)
            thammocks=translate(tk[:hammocks],kcenter)
            #write the posts, then structure then hammocks
            [tposts,tsupport,thammocks]
        end
        #======================old serpentine attempt======================================
        kcols = [kernelmat[i,:] for i in 1:size(kernelmat)[1]]
        #snakify
        for i in 1:length(kcols)
            if iseven(i)
                reverse!(kcols[i])
            end
        end
        #kvec is a vector with structure [[structure1,hammock1],[structure2,hammock2]...]
        kvec = vcat(kcols...)
        ==================================================================================#
        #print the kernels in lexical order. In the future we could switch to a serpentine,
        #but that would require thinking about print direction while generating the kernels
        kvec = reshape(kernelmat,:)
        #return our compiled geometry so we can come back up one level in the directory structure
        #to make the job. This will make the paths nice for our multijob
        #later
        vcat(bumpers,vcat(kvec...))
    end
    GWLJob(joinpath(scaffolddir,"scaffold.gwl"),compgeom...;
           stagespeed=kwargs[:stagespeed],interfacepos=kwargs[:interfacepos])
end

function scaffold(scaffolddir,configfilename::String)
    config = include(configfilename)
    scaffold(scaffolddir,config)
end

scaffold(scaffolddir) = scaffold(scaffolddir,"config.jl")

"""
```julia
arrangescaffolds(arraydims,arrayshape,arraycenter,maxscafdims)
```
Return a `Matrix` with shape `arrayshape` of scaffold coordinate centers centered
on `arraycenter`. Given a maximum scaffold dimension of `maxscafdims` these scaffolds
will be guarenteed not to overlap and to fit in a bounding box with size `arraydims`.
"""
function arrangescaffolds(arraydims,arrayshape,arraycenter,maxscafdims)
    #get the corners of our bounding box
    bboxtopleft = arraycenter + [-arraydims[1], arraydims[2]]/2
    bboxbottomright = arraycenter + [arraydims[1], -arraydims[2]]/2

    #get the coordinates of our top left and bottom right scaffold
    firstcenter = bboxtopleft + [maxscafdims[1],-maxscafdims[2]]/2
    lastcenter = bboxbottomright + [-maxscafdims[1],maxscafdims[2]]/2
    #this is enough to build the matrix
    (xrange,yrange) = map(1:2) do d
        range(start=firstcenter[d],stop=lastcenter[d],length=arrayshape[d])
    end
    @assert step(xrange) > maxscafdims[1] "scaffolds would overlap in x direction"
    @assert (-step(yrange)) > maxscafdims[2] "scaffolds would overlap in y direction"
    centers = [[x,y] for x in xrange, y in yrange]
end

#build a multijob from a matrix of `centercoords => config` pairs
function buildmj(jobs::Matrix{<:Pair})
    #snake it
    stagespeed = nothing
    rowjobs = map(1:size(jobs)[2]) do j
        thisrow = jobs[:,j]
        if iseven(j)
            thisrow = reverse(thisrow)
        end
        thesejobs = map(1:length(thisrow)) do i
            (center,config) = thisrow[i]
            #assume stagespeed is always the same
            stagespeed = config[:stagespeed]
            thisjob = (center => scaffold("$i,$j",config))
            #write the configuration into the scaffold folder so we can look at it later
            open(joinpath("$i,$j","config.txt"), w) do io
                print(io,config)
            end
            return thisjob
        end
    end
    multijob("psweep.gwl",vcat(rowjobs...)...;stagespeed)
end

"""
```julia
psweep([config,] p1 => values[, p2 => values]; arraydims, arrayshape, arraycenter)
```
Build a `multijob` which builds scaffolds with varying parameters. The final array will have
shape `arrayshape` centered on `arraycenter`. These scaffolds will be guarenteed not to overlap
and to fit in a bounding box with size `arraydims`. `config` (provided as a filepath or `Dict`)
should contain all other configuration parameters. If `config` is not provided, a configuration
file is assumed to exist at `"config.jl"`. If swept parameters are included in `config` they will
be ignored.
"""
function psweep end

#for one parameter
function psweep(config::Dict,sweep::Pair{Symbol,<:Vector}; arraydims,arrayshape,arraycenter)
    #destructure our parameter and values
    (p,values) = sweep
    #build a vector of configurations reflecting our parameter sweep
    configs = map(values) do v
        thisconfig = copy(config)
        thisconfig[p] = v
        return thisconfig
    end

    maxscafdims = map([:lscaf, :wscaf]) do dim
        maximum(configs) do c
            c[dim]
        end
    end
    scafcenters = arrangescaffolds(arraydims,arrayshape,arraycenter,maxscafdims)
    @assert length(configs) == length(scafcenters) "Number of parameter values must match array size"
    jobmat = map(zip(scafcenters,reshape(configs,size(scafcenters)...))) do (sc,c)
        sc => c
    end
    buildmj(jobmat)
end

#for two parameters
function psweep(config::Dict,sweep1::Pair{Symbol,<:Vector},sweep2::Pair{Symbol,<:Vector};
                arraydims,arrayshape,arraycenter)
    #destructure our parameter and values
    (p1,values1) = sweep1
    (p2,values2) = sweep2
    #build a matrix representing our parameter combos
    pmat = [Dict(p1 => v1, p2 => v2) for v1 in values1, v2 in values2]
    #now make a matrix of configs
    configs = map(pmat) do pdict
        thisconfig = copy(config)
        for (p,v) in pdict
            thisconfig[p] = v
        end
        return thisconfig
    end

    maxscafdims = map([:lscaf, :wscaf]) do dim
        maximum(configs) do c
            c[dim]
        end
    end
    scafcenters = arrangescaffolds(arraydims,arrayshape,arraycenter,maxscafdims)
    @assert length(configs) == length(scafcenters) "Number of parameter values must match array size"
    jobmat = map(zip(scafcenters,configs)) do (sc,c)
        sc => c
    end
    buildmj(jobmat)
end

function psweep(config::String,args...;kwargs...)
    cdict = include(config)
    psweep(cdict,args...;kwargs...)
end

psweep(args::Vararg{<:Pair{Symbol,<:Vector}};kwargs...) = psweep("config.jl",args...;kwargs...)

end # module Waffle
