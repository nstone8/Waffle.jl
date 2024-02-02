module Waffle
using Tessen, Delica, Unitful
#go get some internal Tessen stuff
import Tessen:HatchLine, pointalong, intersections

function createconfig(filename="config.jl")
    config = Dict(
        :lscaf => 6u"mm",
        :wscaf => 4u"mm",
        :dfield => 1750u"µm",
        :hbottom => 30u"µm",
        :htop => 50u"µm",
        :chamferbottom => pi/6,
        :chamfertop => pi/6,
        :chamferbottom => pi/6,
        :cutangle => pi/6,
        :overlap => 10u"µm",
        :dslice => 1u"µm",
        :wbumper => 200u"µm",
        :fillet => 30u"µm",
        :wpost => 100u"µm",
        :hbeam => 10u"µm",
        :wbeam => 25u"µm",
        :lbeammax => 120u"µm",
        :maxseglength => 50u"µm"
    )
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
#args for current version: endx=(kwargs[:lscaf]/2 - lseg + zmid*tan(kwargs[:cutangle]))
#startx = (-kwargs[:lscaf]/2 + lseg)
function bumpersegtrain(startx,endx,startz,endz,bec;kwargs...)
    #to connect the left and right sides, we will use a bunch of identical segments cut to
    #have overhangs on the left side and placed to overlap with each other
    #get the amount of length that should be devoted to overhangs
    lseg = sqrt(kwargs[:dfield]^2 - kwargs[:wbumper]^2)
    olength = 2*(endz-startz)*tan(kwargs[:cutangle])
    #the length of every slice should be constant, the maximum progress we can make per slice
    #is set by
    lslicemax = lseg - (olength/2) - kwargs[:overlap]
    #the distance we have to cover (including overlapping at the far end is...
    reqdist = endx -  startx + sign(endx-startx)*kwargs[:overlap]

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
    #all the centers of the arc making up the left face will be in the same place
    cleftarc = [-lseg/2 + kwargs[:wbumper]/2, 0u"µm"]
    lbslices = map(range(start=0u"µm",step=kwargs[:dslice],stop=zmid)) do z
        #get the x coordinate of our bottom left corner, this should be -lseg/2
        #for our widest slice
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom right corner, this should be lseg/2
        #for the bottom slice
        xbr = (lseg/2) - z*tan(kwargs[:cutangle])
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
    crightarc = [lseg/2 - kwargs[:wbumper]/2, 0u"µm"]
    rbslices = map(range(start=0u"µm",step=kwargs[:dslice],stop=zmid)) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom left corner, this should be -lseg/2
        #for the topmost slice
        xbl = (-lseg/2) + (zmid-z)*tan(kwargs[:cutangle])
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
    crightarc2 = [lseg/4 - kwargs[:wbumper]/2, 0u"µm"]
    rtslices = map(range(start=zmid2,step=kwargs[:dslice],stop=kwargs[:hbottom]+kwargs[:htop])) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom left corner, this should be -lseg/4
        #for the bottom slice
        xbl = (-lseg/4) + (z)*tan(kwargs[:cutangle])
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

    cleftarc2 = [-lseg/4 + kwargs[:wbumper]/2, 0u"µm"]
    ltslices = map(range(start=zmid2,step=kwargs[:dslice],stop=kwargs[:hbottom]+kwargs[:htop])) do z
        ycoords = bec(z)
        wslice = abs(-(ycoords...))
        #get the x coordinate of our bottom right corner, this should be lseg/4
        #for the topmost slice
        xbr = (lseg/4) + (z-kwargs[:htop])*tan(kwargs[:cutangle])
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
    #the first segment needs to be chamfered according to `chamfertop`, we will do the rest at
    #cutangle
    distperseg = (stopx - startx)/nsegs
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
    #the keystone is the same size as the other segments but cut differently
    keystone = box(lseg,kwargs[:wbeam],kwargs[:hbeam],kwargs[:dslice],
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

#build a kernel with a nx x ny array of posts with pitch of px and py in each dimension
#lbeams and bbeams controls wheter beams should be drawn connecting to geometry to the left
#or below the current kernel
function kernel(nx,ny,px,py,nseg;lbeams=false,bbeams=false,kwargs...)
    #calculate beam lengths
    lbeamx = px - kwargs[:wpost]
    lbeamy = py - kwargs[:wpost]
    #with the optional beams to the left and bottom, the total size of the array is...
    kernelsize = [nx*(kwargs[:wpost] + lbeamx),
                  ny*(kwargs[:wpost] + lbeamy) + lbeamy]
    #if we think of the upper left corner of this maximal array as (0,0) the post coordinates are...
    topleftpostcoords = [[i*(lbeamx + kwargs[:wpost]) - kwargs[:wpost]/2,
                          j*(lbeamy + kwargs[:wpost]) - kwargs[:wpost]/2] for
                             i in 1:nx, j in 1:ny]
    #however, we want the objective at the middle of this maximal array, so let's move all
    #these points to be centered around kernelsize/2
    postcoords = [tlpc - kernelsize/2 for tlpc in topleftpostcoords]
    #little helper function
    function onepost(postcenter,left,bottom)
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
        if left
            push!(allbeams,leftbeam)
        end
        if bottom
            #make the bottom beam by rotating topbeam
            push!(allbeams,rotate(topbeam,pi,postcenter,preserveframe=true))
        end
        #make a namedtuple so we can arrange things nicely in the enclosing scope
        (post=thispost,beams=allbeams)
    end
    needleft = [(i==1) ? lbeams : true for i in 1:nx, j in 1:ny]
    needbottom = [(j==ny) ? bbeams : false for i in 1:nx, j in 1:ny]
    units = map(zip(postcoords,needleft,needbottom)) do (pc,nl,nb)
        onepost(pc,nl,nb)
    end
    uvec = reshape(units,:)
    leftsegs = vcat((b.leftsegs for u in uvec for b in u.beams)...)
    rightsegs = vcat((b.rightsegs for u in uvec for b in u.beams)...)
    keystones = collect(b.keystone for u in uvec for b in u.beams)
    posts = [u.post for u in uvec]
    (posts=posts,leftsegs=leftsegs,rightsegs=rightsegs,keystones=keystones)
end

end # module Waffle
