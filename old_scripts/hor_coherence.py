from collections.abc import Iterable
#from typing import Self
from Bio.Phylo.BaseTree import Tree, Clade
# from Bio.Phylo.PhyloXML import Clade, Phylogeny, Sequence, Phyloxml, Property
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, Location

from loops import LoopInSeq, LoopSpanInSeq
from phylogeny_to_levels import parse_location

def locationContainsLocation(
    containingLocation: SimpleLocation,
    containedLocation: SimpleLocation
) -> bool:
    return (
        containedLocation.start >= containingLocation.start and 
        containedLocation.end <= containingLocation.end
    )

def locationListContainsLocation(
    containingLocationList: list[SimpleLocation],
    containedLocation: SimpleLocation
) -> bool:
    for containingLocation in containingLocationList:
        if locationContainsLocation(
            containingLocation=containingLocation,
            containedLocation=containedLocation
        ):
            return True
    return False

def checkLocationListContainsLocationList(
    containingLocationList: list[SimpleLocation],
    containedLocationList: list[SimpleLocation],
    containementNotFoundTemplate: str = "HOR location {containedLocation} not contained in parent HOR locations {containingLocationList}"
):
    for containedLocation in containedLocationList:
        if not locationListContainsLocation(
            containingLocationList=containingLocationList,
            containedLocation=containedLocation
        ):
            print(containementNotFoundTemplate.format(
                containingLocationList=containingLocationList,
                containedLocation=containedLocation
            ))
            
def locationOverlapsLocation(
    location1: SimpleLocation,
    location2: SimpleLocation
) -> bool:
    return (
        location1.start < location2.end and 
        location1.end > location2.start
    )

def locationListOverlapsLocationList(
    locationList1: list[SimpleLocation],
    locationList2: list[SimpleLocation]
) -> bool:
    for location1 in locationList1:
        for location2 in locationList2:
            if locationOverlapsLocation(
                location1=location1,
                location2=location2
            ):
                return True
    return False

def checkLocationListOverlapsLocationList(
    locationList1: list[SimpleLocation],
    locationList2: list[SimpleLocation],
    foundOverlapTemplate: str = "Found overlap between {location1} and {location2}"
):
    for location1 in locationList1:
        for location2 in locationList2:
            if locationOverlapsLocation(
                location1=location1,
                location2=location2
            ):
                print(foundOverlapTemplate.format(
                    location1=location1,
                    location2=location2
                ))

def checkLocationListSelfOverlap(
    locationList: list[SimpleLocation],
    foundOverlapTemplate: str = "Found overlap between {location1} and {location2}"
):
    controlledLocationList = []
    for location in locationList:
        checkLocationListOverlapsLocationList(
            locationList1=controlledLocationList,
            locationList2=[location],
            foundOverlapTemplate=foundOverlapTemplate
        )
        controlledLocationList.append(location)

def checkHORCoherence(horClade: Clade, path: list[int]):
    parentLocations: list[SimpleLocation] = [parse_location(seq.location) for seq in horClade.sequences]
    horName = horClade.name if horClade.name is not None else 'unknown'
    checkLocationListSelfOverlap(
        parentLocations, 
        foundOverlapTemplate=(
            f"In clade {horName} (path: {path}) " +
            "found self overlap between {location1} and {location2}"
        )
    )
    locationsFound: list[SimpleLocation] = []
    for subhor in horClade.clades:
        subhorName = subhor.name if subhor.name is not None else 'unknown'
        newLocations = [parse_location(seq.location) for seq in subhor.sequences]
        checkLocationListContainsLocationList(
            containingLocationList=parentLocations,
            containedLocationList=newLocations,
            containementNotFoundTemplate=(
                f"Location of HOR {subhorName} " +
                "({containedLocation}) " +
                f"not found among locations of HOR {horName} " +
                "({containingLocationList})"
            )
        )
        checkLocationListOverlapsLocationList(
            locationList1=locationsFound,
            locationList2=newLocations,
            foundOverlapTemplate=(
                f"In clade {subhorName} " +
                "found sibling overlap between {location1} and {location2}"
            )
        )
        locationsFound.extend(newLocations)
    for subhorIndex, subhor in enumerate(horClade.clades):
        checkHORCoherence(subhor, path=path + [subhorIndex])
        
def checkHORTreeCoherence(horTree: Tree):
    checkHORCoherence(horClade=horTree.root, path=[])
    
    
def loopSpanInSeqOverlapsLoopSpanInSeq(
    loopSpanInSeq1: LoopSpanInSeq,
    loopSpanInSeq2: LoopSpanInSeq
) -> bool:
    return (
        loopSpanInSeq1.span_start < loopSpanInSeq2.span_start + loopSpanInSeq2.span_length and 
        loopSpanInSeq1.span_start + loopSpanInSeq1.span_length > loopSpanInSeq2.span_start
    )
         
def checkSpanInSeqListOverlapsSpanInSeqList(
    loopSpanInSeqList1: list[LoopSpanInSeq],
    loopSpanInSeqList2: list[LoopSpanInSeq],
    foundOverlapTemplate: str = "Found overlap between {loopSpanInSeq1} and {loopSpanInSeq2}"
):
    for loopSpanInSeq1 in loopSpanInSeqList1:
        for loopSpanInSeq2 in loopSpanInSeqList2:
            if loopSpanInSeqOverlapsLoopSpanInSeq(
                loopSpanInSeq1=loopSpanInSeq1,
                loopSpanInSeq2=loopSpanInSeq2
            ):
                print(foundOverlapTemplate.format(
                    loopSpanInSeq1=loopSpanInSeq1,
                    loopSpanInSeq2=loopSpanInSeq2
                ))

def checkLoopInSeqSelfOverlap(
    loopInSeq: LoopInSeq,
    foundOverlapTemplate: str = "Found overlap between {loopSpanInSeq1} and {loopSpanInSeq2}"
):
    controlledSpanInSeqList = []
    for spanInSeq in loopInSeq.spans_in_seq:
        checkSpanInSeqListOverlapsSpanInSeqList(
            loopSpanInSeqList1=controlledSpanInSeqList,
            loopSpanInSeqList2=[spanInSeq],
            foundOverlapTemplate=foundOverlapTemplate
        )
        controlledSpanInSeqList.append(spanInSeq)
            
# def checkLoopInSeqOverlapsLoopInSeq(
#     loopInSeq1: LoopInSeq,
#     loopInSeq2: LoopInSeq
# ) -> bool:
#     return (
#         location1.start < location2.end and 
#         location1.end > location2.start
#     )


def spanOperlapsSpan(span1: tuple[int,int], span2: tuple[int,int]):
    return (
        span1[0] < span2[1] and 
        span1[1] > span2[0]
    )

def checkSpanListOverlapsSpanList(
    spanList1: list[tuple[int,int]],
    spanList2: list[tuple[int,int]],
    foundOverlapTemplate: str = "Found overlap between {span1} and {span2}"
):
    for span1 in spanList1:
        for span2 in spanList2:
            if spanOperlapsSpan(
                span1=span1,
                span2=span2
            ):
                print(foundOverlapTemplate.format(
                    span1=span1,
                    span2=span2
                ))

def checkSpanListSelfOverlap(
    spanList: list[tuple[int,int]],
    foundOverlapTemplate: str = "Found overlap between {loopSpanInSeq1} and {loopSpanInSeq2}"
):
    controlledSpanList = []
    for span in spanList:
        checkSpanListOverlapsSpanList(
            spanList1=controlledSpanList,
            spanList2=[span],
            foundOverlapTemplate=foundOverlapTemplate
        )
        controlledSpanList.append(span)
