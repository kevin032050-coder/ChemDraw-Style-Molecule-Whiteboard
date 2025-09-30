from cmu_graphics import *
import requests
import math
import copy
from rdkit import Chem
from rdkit.Chem import Draw
import os
from rdkit.Chem.EnumerateStereoisomers import *
from rdkit.Chem import rdMolDescriptors

def distance(x1, y1, x2, y2):
    return ((x1-x2)**2+(y1-y2)**2)**0.5

def createElement(s):
    possibleElements = ['Carbon', 'Oxygen', 'Nitrogen', 'Hydrogen', 'Phosphorus' 
                        ,'Sulfur', 'Fluorine', 'Chlorine', 'Bromine', 'Iodine']
    if s == 'Carbon':
        return Carbon(None, None)
    elif s == 'Oxygen':
        return Oxygen(None, None)
    elif s == 'Nitrogen':
        return Nitrogen(None, None)
    elif s=='Hydrogen':
        return Hydrogen(None, None)
    elif s=='Phosphorus':
        return Phosphorus(None, None)
    elif s=='Sulfur':
        return Sulfur(None, None)
    elif s=='Fluorine':
        return Fluorine(None, None)
    elif s=='Chlorine':
        return Chlorine(None, None)
    elif s=='Bromine':
        return Bromine(None, None)
    elif s=='Iodine':
        return Iodine(None, None)

def generateImage(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Draw.MolToFile(mol, "image.png")
    return 'image.png'

def generateSmilesHelper(app, start_element):
    # reset visited status 
    for element in app.selectedElements:
        element.visited = False
        element.ringIndex = None  

    ringTracker = {}  
    smiles = search(start_element, None, ringTracker)
    return smiles

def search(element, parent, ringTracker):
    element.visited = True
    smiles = element.name[0].upper()  

    branches = []
    for bond in element.bonds:
        neighbor = bond.element2 if bond.element1 == element else bond.element1

        if neighbor == parent:
            continue  # skip the bond to the parent element

        bondChar = "" 
        if bond.bondType == "double":
            bondChar = "="
        elif bond.bondType == "triple":
            bondChar = "#"

        if not neighbor.visited:
            branches.append(bondChar + search(neighbor, element, ringTracker))
        else:
            # handle ring closure
            if neighbor.ringIndex is None:
                # assign a new ring closure index
                ringIndex = len(ringTracker) + 1
                neighbor.ringIndex = ringIndex
                element.ringIndex = ringIndex
                ringTracker[ringIndex] = True
                smiles += f"{bondChar}{ringIndex}"
            else:
                # use the existing ring index
                smiles += f"{bondChar}{neighbor.ringIndex}"

    # add branches to SMILES
    if len(branches) > 1:
        smiles += "".join(f"({branch})" for branch in branches[:-1])
        smiles += branches[-1]
    elif len(branches) == 1:
        smiles += branches[0]

    return smiles

def findStereocenters(elements, bonds):
    stereocenters = {"chiralCenters": [], "doubleBonds": [] }

    for element in elements:
        if len(element.bonds) == 4:
            substituents = [bond.element2 if bond.element1 == element else 
                            bond.element1 for bond in element.bonds]
            if len(set(substituents)) == 4:
                stereocenters["chiralCenters"].append(element)
    for bond in bonds:
        if bond.bondType == "double":
            neighbors1 = []
            for b in bond.element1.bonds:
                if b != bond:
                    if b.element1 == bond.element1:
                        neighbors1.append(b.element2)
                    else:
                        neighbors1.append(b.element1)
            neighbors2 = []
            for b in bond.element2.bonds:
                if b != bond:
                    if b.element1 == bond.element2:
                        neighbors2.append(b.element2)
                    else:
                        neighbors2.append(b.element1)
            if len(set(neighbors1)) > 0 and len(set(neighbors2)) > 0:
                stereocenters["doubleBonds"].append(bond)
    return stereocenters

def enumerateStereoisomers(elements, bonds, stereocenters, index=0, 
                            currentConfig=None, allConfigs=None):

    if currentConfig is None:
        currentConfig = []
    if allConfigs is None:
        allConfigs = []
    #base case
    if index == (len(stereocenters["chiralCenters"]) + 
                len(stereocenters["doubleBonds"])):
        allConfigs.append(copy.copy(currentConfig))
        return allConfigs

    #recursive case
    if index < len(stereocenters["chiralCenters"]):
        #backtracking
        for stereo in ["@", "@@"]:
            currentConfig.append(stereo)
            enumerateStereoisomers(elements, bonds, stereocenters, index + 1, 
                                    currentConfig, allConfigs)
            currentConfig.pop()
    else:
        for stereo in ["/", "\\"]:
            currentConfig.append(stereo)
            enumerateStereoisomers(elements, bonds, stereocenters, index + 1, 
                                    currentConfig, allConfigs)
            currentConfig.pop()

    return allConfigs

def applyStereochemistry(elements, bonds, stereocenters, configuration):
    for i in range(len(stereocenters["chiralCenters"])):
        center = stereocenters["chiralCenters"][i]
        stereo = configuration[i]
        center.stereo = stereo

    for i in range(len(stereocenters["doubleBonds"])):
        bond = stereocenters["doubleBonds"][i]
        stereo = configuration[len(stereocenters["chiralCenters"]) + i]
        bond.stereo = stereo

def resetStereochemistry(elements, bonds):
    for element in elements:
        element.stereo = None
    for bond in bonds:
        bond.stereo = None

def generateSmilesWithStereo(elements, startElement, parent=None, visited=None):
    if visited is None:
        visited = set()

    visited.add(startElement)
    smiles = startElement.name[0].upper()
    branches = []
    for bond in startElement.bonds:
        neighbor = (bond.element2 if bond.element1 == startElement 
                    else bond.element1)

        if neighbor == parent:
            continue

        if neighbor in visited:
            ringIdx = bond.ringIndex if hasattr(bond, "ringIndex") else ""
            bondChar = "" if bond.bondType == "single" else bond.bondType[0]
            branches.append(bondChar + str(ringIdx))
            continue
        bondChar = ""
        stereoPrefix = ""
        if bond.bondType == "double":
            bondChar = "="
            if bond.stereo == "/":
                stereoPrefix = "/"
            elif bond.stereo == "\\":
                stereoPrefix = "\x5c"
        branchSmiles = generateSmilesWithStereo(elements, neighbor, startElement, 
                                                visited)
        if bond.bondType == "double" and stereoPrefix:
            branchSmiles = f"{stereoPrefix}{branchSmiles}"
        branches.append(bondChar + branchSmiles)

    if len(branches) > 1:
        smiles += ''.join(f"({branch})" for branch in branches[:-1])
        smiles += branches[-1]
    elif branches:
        smiles += branches[0]

    return smiles

def smilesWithStereochemistry(app):
    # used chatgpt with logic help
    # find stereocenters
    stereocenters = findStereocenters(app.selectedElements, app.bonds)

    # generate stereoisomer configurations
    configurations = enumerateStereoisomers(app.selectedElements, app.bonds, 
                                            stereocenters)

    # generate SMILES for each configuration
    stereoisomers = []
    for config in configurations:
        resetStereochemistry(app.selectedElements, app.bonds)  
        applyStereochemistry(app.selectedElements, app.bonds, stereocenters, 
                            config)  
        smiles = generateSmilesWithStereo(app.selectedElements, 
                                        app.selectedElements[0])  
        stereoisomers.append(smiles)

    return stereoisomers

def getIUPACName(smiles):
    # PubChem API
    url = url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/" +
        "property/IUPACName/TXT")
    response = requests.get(url)
    
    if response.status_code == 200:
        name = response.text.strip()  # he IUPAC name
        index = 0
        while index < len(name) and not name[index].isalpha():
            index += 1
        if index < len(name):  
            name = name[:index] + name[index].upper() + name[index + 1:]
        return name
    else:
        return None

def getPropertyInfo(smiles):
    properties = ['IUPACName', 'MolecularFormula', 'MolecularWeight', 
    'HBondDonorCount', 'HBondAcceptorCount']
    result = []
    for property in properties:
        url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/" +
        f"property/{property}/TXT")
        response = requests.get(url)
        if response.status_code == 200:
            if property == 'IUPACName':
                name = response.text.strip()  # he IUPAC name
                index = 0
                while index < len(name) and not name[index].isalpha():
                    index += 1
                if index < len(name):  
                    name = name[:index] + name[index].upper() + name[index + 1:]
                result.append(name)
            else:
                result.append(response.text.strip())
        else:
            result.append(None)
    return result

class Element:
    def __init__(self):
        self.visited = False
        self.ringIndex = None
        self.stereo = None
    def addBond(self, bondType):
        bond = Bond(bondType, self)
        self.bonds.append(bond)
        return bond
    
    def __repr__(self):
        return f'name={self.name}, bonds={self.bonds}'
    
    def findBondNum(self):
        count = 0
        for bond in self.bonds:
            count += bond.getBondCount()
        return count

class Carbon(Element):
    carbonNum = 0
    def __init__(self, x, y):
        Carbon.carbonNum += 1
        self.name = 'carbon'+str(Carbon.carbonNum)
        self.bonds = []
        self.maxBonds = 4
        self.x = x
        self.y = y

class Nitrogen(Element):
    nitrogenNum = 0
    def __init__(self, x, y):
        Nitrogen.nitrogenNum += 1
        self.name = 'nitrogen'+str(Nitrogen.nitrogenNum)
        self.bonds = []
        self.maxBonds = 3
        self.x = x
        self.y = y

class Oxygen(Element):
    oxygenNum = 0
    def __init__(self, x, y):
        Oxygen.oxygenNum += 1
        self.name = 'oxygen'+str(Oxygen.oxygenNum)
        self.bonds = []
        self.maxBonds = 2
        self.x = x
        self.y = y

class Hydrogen(Element):
    hydrogenNum = 0
    def __init__(self, x, y):
        Hydrogen.hydrogenNum += 1
        self.name = 'hydrogen'+str(Hydrogen.hydrogenNum)
        self.bonds = []
        self.maxBonds = 1
        self.x = x
        self.y = y

class Phosphorus(Element):
     phosphorusNum = 0
     def __init__(self, x, y):
        Phosphorus.phosphorusNum += 1
        self.name = 'phosphorus'+str(Phosphorus.phosphorusNum)
        self.bonds = []
        self.maxBonds = 4
        self.x = x
        self.y = y

class Sulfur(Element):
    sulfurNum = 0
    def __init__(self, x, y):
        Sulfur.sulfurNum += 1
        self.name = 'sulfur'+str(Sulfur.sulfurNum)
        self.bonds = []
        self.maxBonds = 4
        self.x = x
        self.y = y

class Fluorine(Element):
    fluorineNum = 0
    def __init__(self, x, y):
        Fluorine.fluorineNum += 1
        self.name = 'fluorine'+str(Fluorine.fluorineNum)
        self.bonds = []
        self.maxBonds = 1
        self.x = x
        self.y = y

class Chlorine(Element):
    chlorineNum = 0
    def __init__(self, x, y):
        Chlorine.chlorineNum += 1
        self.name = 'chlorine'+str(Chlorine.chlorineNum)
        self.bonds = []
        self.maxBonds = 1
        self.x = x
        self.y = y

class Bromine(Element):
    bromineNum = 0
    def __init__(self, x, y):
        Bromine.bromineNum += 1
        self.name = 'bromine'+str(Bromine.bromineNum)
        self.bonds = []
        self.maxBonds = 1
        self.x = x
        self.y = y

class Iodine(Element):
    iodineNum = 0
    def __init__(self, x, y):
        Iodine.iodineNum += 1
        self.name = 'iodine'+str(Iodine.iodineNum)
        self.bonds = []
        self.maxBonds = 1
        self.x = x
        self.y = y

class Bond:
    def __init__(self, bondType, element1):
        self.bondType = bondType
        self.element1 = element1
        self.element2 = None
        self.stereo = None

    def addElement(self, element):
        self.element2 = element
        element.bonds.append(self)

    def getBondCount(self):
        if self.bondType == 'single':
            return 1
        elif self.bondType == 'double':
            return 2
        elif self.bondType == 'triple':
            return 3
        return 1
    
    def __repr__(self):
        return (f'name={self.bondType}, element1={self.element1.name},' + 
        f'element2={self.element2.name if self.element2 != None else None}')

class Board():
    def __init__(self, rows, cols, boardLeft, boardTop, boardWidth, boardHeight):
        self.rows = rows
        self.cols = cols
        self.boardLeft = boardLeft
        self.boardTop = boardTop
        self.boardWidth = boardWidth
        self.boardHeight = boardHeight
        self.cellBorderWidth = 2
    
    def drawBoard(self):
        for row in range(self.rows):
            for col in range(self.cols):
                self.drawCell(row, col)
        self.drawBoardBorder()

    def drawBoardBorder(self):
        drawRect(self.boardLeft, self.boardTop, self.boardWidth, self.boardHeight,
            fill=None, border=rgb(38, 84, 124),
            borderWidth=2*self.cellBorderWidth)

    def drawCell(self, row, col):
        cellLeft, cellTop = self.getCellLeftTop(row, col)
        cellWidth, cellHeight = self.getCellSize()
        drawRect(cellLeft, cellTop, cellWidth, cellHeight,
                fill=None, border=rgb(38, 84, 124),
                borderWidth=self.cellBorderWidth)

    def getCellLeftTop(self, row, col):
        cellWidth, cellHeight = self.getCellSize()
        cellLeft = self.boardLeft + col * cellWidth
        cellTop = self.boardTop + row * cellHeight
        return (cellLeft, cellTop)

    def getCellSize(self):
        cellWidth = self.boardWidth / self.cols
        cellHeight = self.boardHeight / self.rows
        return (cellWidth, cellHeight)

def onAppStart(app):
    #canvas size
    app.width = 800
    app.height = 600
    app.SMILES = None
    app.bondSelected = 'single'
    app.maxBondReached = False
    app.stereo = None
    app.ringSelected = None
    app.background = rgb(201, 228, 231)
    app.title = 'Organic Molecule Visualizer'

def start_redrawAll(app):
    drawLabel('Organic Molecule Visualizer', app.width//2, app.height//2, size=60, 
    font='arial')
    drawLabel("Press 'space' to continue", app.width//2, app.height//2+50, size=25, 
    font = 'monospace')

def start_onKeyPress(app, key):
    if key == 'space':
        setActiveScreen('modeSelection')

def modeSelection_onScreenActivate(app):
    app.selected = False
    app.selectedElements = []
    app.points = []
    app.lastElement = None
    app.nextElement = None
    app.bonds = []
    app.SMILES = None

def modeSelection_redrawAll(app):
    drawRect(200, 200, 400, 75, fill=rgb(134, 176, 188), borderWidth=5, 
            border=rgb(38, 84, 124))
    drawRect(200, 300, 400, 75, fill=rgb(134, 176, 188), borderWidth=5, 
            border=rgb(38, 84, 124))
    drawRect(200, 400, 400, 75, fill=rgb(134, 176, 188), borderWidth=5, 
            border=rgb(38, 84, 124))
    drawLabel('Organic Molecule Visualizer', app.width//2, 150, size=40, 
            font='arial')
    drawLabel('Draw Mode', 400, 238, size=30)
    drawLabel('Name Mode', 400, 338, size=30)
    drawLabel('Common Organic Molecules', 400, 438, size=30)

def modeSelection_onMousePress(app, mouseX, mouseY):
    if 200<=mouseX<=600 and 200<=mouseY<=275:
        setActiveScreen('drawPage')
    elif 200<=mouseX<=600 and 300<=mouseY<=375:
        setActiveScreen('namePage')
    elif 200<=mouseX<=600 and 400<=mouseY<=475:
        setActiveScreen('commonPage') 

def drawPage_onScreenActivate(app):
    app.selected = False
    app.selectedElements = []
    app.points = []
    app.lastElement = None
    app.nextElement = None
    app.bonds = []
    app.organicElements = ['Carbon', 'Oxygen', 'Nitrogen', 'Hydrogen', 
                           'Phosphorus', 'Sulfur', 'Fluorine', 'Chlorine', 
                           'Bromine', 'Iodine']
    app.nonLinear = False
    app.nonLinearElement = None
    app.width = 800
    app.height = 600
    app.bondSelected = 'single'
    app.possibleBonds = ['single', 'double', 'triple']
    app.maxBondReached = False
    app.ringSelected = None

def drawPage_redrawAll(app):
    drawLine(75, 0, 75, 600, fill=rgb(38, 84, 124))

    elements = ['C', 'O', 'N', 'H', 'P', 'S', 'F', 'Cl', 'Br', 'I']
    for y in range(10, 556, 60):
        drawRect(15, y, 45, 45, fill=None, border=rgb(38, 84, 124))
        elementNum = (y-10)//60
        drawLabel(f'{elements[elementNum]}', 37.5, y+23, size=25)
    
    drawLine(725, 0, 725, 600, fill=rgb(38, 84, 124))
    
    if app.bondSelected == 'single' or app.bondSelected == None:
        drawRect(740, 15, 45, 45, fill='gray', border=rgb(38, 84, 124), 
                opacity=35)
    elif app.bondSelected == 'double':
        drawRect(740, 75, 45, 45, fill='gray', border=rgb(38, 84, 124), 
                opacity=35)
    elif app.bondSelected == 'triple':
         drawRect(740, 135, 45, 45, fill='gray', border=rgb(38, 84, 124), 
                opacity=35)
    #single bond
    drawRect(740, 15, 45, 45, fill=None, border=rgb(38, 84, 124))
    drawLine(747, 37.5, 778, 37.5)
    #double bond
    drawRect(740, 75, 45, 45, fill=None, border=rgb(38, 84, 124))
    drawLine(747, 95, 778, 95)
    drawLine(747, 100, 778, 100)
    #triple bond
    drawRect(740, 135, 45, 45, fill=None, border=rgb(38, 84, 124))
    drawLine(747, 152.5, 778, 152.5)
    drawLine(747, 157.5, 778, 157.5)
    drawLine(747, 162.5, 778, 162.5)


    for y in range(195, 371, 60):
        drawRect(740, y, 45, 45, fill=None, border=rgb(38, 84, 124))
        points = (y-195)//60
        if points==2:
            drawRegularPolygon(763, y+22.5, 18, 6-points, fill=None, 
                                border='black', rotateAngle=45)
        else:
            drawRegularPolygon(763, y+22.5, 18, 6-points, fill=None, border='black')

    drawRect(740, 550, 45, 45, fill=None, border=rgb(38, 84, 124))
    drawLabel('next', 763, 573, size=15)
    
    if len(app.selectedElements)==1:
        x0, y0 = app.selectedElements[0].x, app.selectedElements[0].y
        if x0 != None:
            drawLabel(f'{app.selectedElements[0].name[0].upper()}', x0, y0, size=15)
    elif len(app.selectedElements)==2 and app.selectedElements[1].x == None:
        x0, y0 = app.selectedElements[0].x, app.selectedElements[0].y
        drawLabel(f'{app.selectedElements[0].name[0].upper()}', x0, y0, size=15)
    else:
        for bond in app.bonds:
            element1 = bond.element1
            element2 = bond.element2
            if element1.x != None:
                x0, y0 = element1.x, element1.y
            else:
                x0 = None
            if element2 != None:
                if element2.x != None:
                    x1, y1 = element2.x, element2.y
                else:
                    x1 = None
            else:
                x1 = None
            if x1 != None:
                drawBond(x0, y0, x1, y1, bond.bondType)
                drawLabel(f'{element1.name[0].upper()}', x0, y0, size=15)
                drawLabel(f'{element2.name[0].upper()}', x1, y1, size=15)
            elif x1==None and x0 != None:
                drawCircle(x0, y0, 7, fill='black')
                drawLabel(f'{element1.name}', x0, y0-12, size=15)
            
    if app.nonLinearElement != None:
        element = app.selectedElements[app.nonLinearElement]
        x, y = element.x, element.y
        drawCircle(x, y, 10, fill=None, border='red')
    elif len(app.selectedElements)>0:
        if app.nextElement == None:
            if app.lastElement != None:
                element = app.lastElement
            else:
                element = None
        else:
            element = app.nextElement
        if element != None:
            x, y = element.x, element.y
            drawCircle(x, y, 10, fill=None, border='red')
    
    if app.maxBondReached == True:
        drawRect(0, 0, 800, 800, fill='gray', opacity=35)
        drawLabel('Max Bonds Reached', app.width//2, app.height//2, size=25)
        drawLabel('Select New Element to Continue', app.width//2, 
                    app.height//2+30, size=23)

def drawBond(x0, y0, x1, y1, bondType):
    #direction vector
    dx = x1 - x0
    dy = y1 - y0
    length = math.sqrt(dx**2 + dy**2)
    if length == 0:
        return  

    # normalize
    dx /= length
    dy /= length

    # adjust the start and end points
    distanceOffset = 6
    startX = x0 + dx * distanceOffset
    startY = y0 + dy * distanceOffset
    endX = x1 - dx * distanceOffset
    endY = y1 - dy * distanceOffset

    if bondType == 'single':
        drawLine(startX, startY, endX, endY, lineWidth=1)
    elif bondType == 'double' or bondType == 'triple':
        # Perpendicular offset vector
        offsetX = -dy
        offsetY = dx

        # Spacing between parallel lines
        spacing = 3
        offsetX *= spacing
        offsetY *= spacing

        if bondType == 'double':
            # Two parallel lines for double bond
            drawLine(startX + offsetX, startY + offsetY, endX + offsetX, 
                    endY + offsetY, lineWidth=1)
            drawLine(startX - offsetX, startY - offsetY, endX - offsetX, 
                    endY - offsetY, lineWidth=1)
        elif bondType == 'triple':
            # Central line for triple bond
            drawLine(startX, startY, endX, endY, lineWidth=1)
            # Offset lines
            drawLine(startX + offsetX, startY + offsetY, endX + offsetX, 
                    endY + offsetY, lineWidth=1)
            drawLine(startX - offsetX, startY - offsetY, endX - offsetX, 
                    endY - offsetY, lineWidth=1)  

def drawPage_onMousePress(app, mouseX, mouseY):
    #max bond logic
    if app.maxBondReached:
        app.maxBondReached = False

    # placing a ring
    if app.selected == True and app.ringSelected != None and 75<mouseX<700:
        placeNewRing(app, mouseX, mouseY)

    # placing new element
    elif app.selected == True and 75<mouseX<725:
        placeNewElement(app, mouseX, mouseY)

    
    # check if selecting new element
    checkSelectingNewElement(app, mouseX, mouseY)
    
    #check if making a ring structure
    if app.selected == False and app.bondSelected != None:
        checkMakingRingStructure(app, mouseX, mouseY)
        
    # check if selecting already drawn elements
    if app.selected == False and mouseX>75:
        checkSelectingDrawnElement(app, mouseX, mouseY)
        
    #check if selecting bonds
    checkSelectingBond(app, mouseX, mouseY)

    #check if selecting ring structure
    checkSelectingRingStructure(app, mouseX, mouseY)

    checkGoNext(app, mouseX, mouseY)
            
def placeNewElement(app, mouseX, mouseY):
    #first element
    if len(app.selectedElements)==1:
        app.lastElement = app.selectedElements[-1]
        x1, y1 = mouseX, mouseY
        app.selectedElements[-1].x = x1
        app.selectedElements[-1].y = y1
        app.points.append((x1, y1))
    #second element
    elif len(app.selectedElements)==2:
        app.nextElement = app.selectedElements[-1]
        if app.bondSelected == None:
            app.bondSelected = 'single'
        app.lastElement.addBond(app.bondSelected)
        app.bonds.append(app.lastElement.bonds[-1])
        app.lastElement.bonds[-1].addElement(app.nextElement)
        x0, y0 = app.lastElement.x, app.lastElement.y
        x1, y1 = calculateNewPosition(x0, y0, mouseX, mouseY, 50)
        app.selectedElements[-1].x = x1
        app.selectedElements[-1].y = y1 
        app.points.append((x1, y1))
    #rest
    elif len(app.selectedElements)>2:
        if app.nonLinear == False:
            app.lastElement = app.nextElement
            app.nextElement = app.selectedElements[-1]
            if app.bondSelected == None:
                app.bondSelected = 'single'
            app.lastElement.addBond(app.bondSelected)
            if app.lastElement.findBondNum() > app.lastElement.maxBonds:
                app.lastElement.bonds.pop()
                app.selectedElements.pop()
                app.nextElement = None
                app.maxBondReached = True
            else:
                app.bonds.append(app.lastElement.bonds[-1])
                app.lastElement.bonds[-1].addElement(app.nextElement)
                x0, y0 = app.lastElement.x, app.lastElement.y
                x1, y1 = calculateNewPosition(x0, y0, mouseX, mouseY, 50)
                app.selectedElements[-1].x = x1
                app.selectedElements[-1].y = y1 
                app.points.append((x1, y1))
        else:
            app.nextElement = app.selectedElements[-1] 
            if app.bondSelected == None:
                app.bondSelected = 'single'
            app.lastElement.addBond(app.bondSelected)
            if app.lastElement.findBondNum() > app.lastElement.maxBonds:
                app.lastElement.bonds.pop()
                app.selectedElements.pop()
                app.nextElement = None
                app.maxBondReached = True
            else:
                app.bonds.append(app.lastElement.bonds[-1])
                app.lastElement.bonds[-1].addElement(app.nextElement)
                x0, y0 = app.lastElement.x, app.lastElement.y
                x1, y1 = calculateNewPosition(x0, y0, mouseX, mouseY, 50)
                app.selectedElements[-1].x = x1
                app.selectedElements[-1].y = y1 
                app.points.append((x1, y1))
                app.nonLinear = False


    app.selected = False
    app.nonLinearElement = None
    app.bondSelected = None
    return

def placeNewRing(app, mouseX, mouseY):
    if app.ringSelected is not None:
        ringSize = app.ringSelected  
        centerX, centerY = mouseX, mouseY
        radius = 40  
        angleStep = 360 / ringSize
        #create the ring elements centered at mouseX, mouseY
        initialRingPoints = []
        for i in range(ringSize):
            angle = math.radians(i * angleStep)
            x = centerX + radius * math.cos(angle)
            y = centerY + radius * math.sin(angle)
            initialRingPoints.append((x, y))

        # if lastElement exists
        closestPointIndex = 0  
        if app.nextElement != None:
            app.lastElement = app.nextElement
        if app.lastElement is not None:
            x0, y0 = app.lastElement.x, app.lastElement.y
            closestDistance = float("inf")
            for i in range(ringSize):
                x, y = initialRingPoints[i]
                dist = distance(x, y, x0, y0)
                if dist < closestDistance:
                    closestDistance = dist
                    closestPointIndex = i
            closestX, closestY = initialRingPoints[closestPointIndex]

            # calculate the new position 
            adjustedClosestX, adjustedClosestY = calculateNewPosition(x0, y0, 
                                                closestX, closestY, 50)
            # calculate shift
            deltaX = adjustedClosestX - closestX
            deltaY = adjustedClosestY - closestY
            # shift all points 
            adjustedRingPoints = [(x + deltaX, y + deltaY) for x, y in 
                                    initialRingPoints]
        else:
            # no last element
            adjustedRingPoints = initialRingPoints

        # create elements and set positions
        ringElements = []
        for x, y in adjustedRingPoints:
            element = createElement('Carbon')  
            element.x, element.y = x, y
            ringElements.append(element)
            app.selectedElements.append(element)
            app.points.append((x, y))

        # connect the last element to the closest ring element
        if app.lastElement is not None:
            bond = app.lastElement.addBond('single')
            bond.addElement(ringElements[closestPointIndex])
            app.bonds.append(bond)

        # create bonds between elements in the ring
        for i in range(ringSize):
            element1 = ringElements[i]
            element2 = ringElements[(i + 1) % ringSize] 
            bond = element1.addBond('single')
            bond.addElement(element2)
            app.bonds.append(bond)

        # find the furthest element from lastElement 
        if app.lastElement is not None:
            x0, y0 = app.lastElement.x, app.lastElement.y
            furthestDistance = float("-inf")
            furthestPointIndex = 0
            for i in range(ringSize):
                x, y = adjustedRingPoints[i]
                dist = distance(x, y, x0, y0)
                if dist > furthestDistance:
                    furthestDistance = dist
                    furthestPointIndex = i
            app.nextElement = ringElements[furthestPointIndex]
        else:
            app.nextElement = ringElements[-1]  
        
        app.ringSelected = None 
        return

def calculateNewPosition(x0, y0, mouseX, mouseY, distance):
    dx = mouseX - x0
    dy = mouseY - y0
    length = (dx**2 + dy**2)**0.5
    if length == 0:
        return x0, y0
    dx = (dx / length) * distance
    dy = (dy / length) * distance
    newX = x0 + dx
    newY = y0 + dy
    return newX, newY

def checkSelectingNewElement(app, mouseX, mouseY):
    for y in range(10, 556, 60):
        y0, y1 = y, y + 45
        newElement = app.organicElements[(y - 10) // 60]
        if 15 <= mouseX <= 60 and y0 <= mouseY <= y1:
            app.selected = True
            app.selectedElements.append(createElement(newElement))
            return

def checkSelectingDrawnElement(app, mouseX, mouseY):
    if len(app.points)==1:
        return 
    for i in range(len(app.points)):
        x, y = app.points[i]
        if distance(x, y, mouseX, mouseY)<7:
            app.lastElement = app.selectedElements[i]
            app.nonLinear = True
            app.nonLinearElement = i
            return

def checkSelectingBond(app, mouseX, mouseY):
    for y in range(15, 181, 60):
        y0, y1 = y, y+45
        newBond = app.possibleBonds[(y-15)//60]
        if 740<=mouseX<=785 and y0 <=mouseY<=y1:
            app.bondSelected = newBond
            return

def checkMakingRingStructure(app, mouseX, mouseY):
    if len(app.points)==1:
            return
    for i in range(len(app.points)):
        x, y = app.points[i]
        if distance(x, y, mouseX, mouseY)<7:
            app.lastElement = app.nextElement
            app.nextElement = app.selectedElements[i]
            app.lastElement.addBond(app.bondSelected)
            app.lastElement.bonds[-1].addElement(app.nextElement)
            if( (app.lastElement.findBondNum() > app.lastElement.maxBonds) or 
            (app.nextElement.findBondNum() > app.nextElement.maxBonds) ):
                app.lastElement.bonds.pop()
                app.nextElement.bonds.pop()
                app.maxBondReached = True
            else:
                app.bonds.append(app.lastElement.bonds[-1])
                app.bondSelected = None

def checkSelectingRingStructure(app, mouseX, mouseY):
    for y in range(195, 371, 60):
        points = 6-(y-195)//60
        y0, y1 = y, y+45
        if 740<=mouseX<=785 and y0 <=mouseY<=y1:
            app.ringSelected = points
            app.selected = True
            return

def checkGoNext(app, mouseX, mouseY):
    if 740<mouseX<740+45 and 550<mouseY<550+45:
        setActiveScreen('resultPage')

def drawPage_onKeyPress(app, key):
    if key == 's':
        generateSmiles(app)
    elif key == 'space':
        app.SMILES = generateSmiles(app)
        setActiveScreen('resultPage')

def namePage_onScreenActivate(app):
    app.name = ''

def namePage_redrawAll(app):
    drawLabel('Enter IUPAC Name:', app.width//2, 200, size=27, bold=True)
    if app.name=='':
        drawLabel('Type to enter name', app.width//2, 250, size=23, italic=True)
    else:
        drawLabel(f'{app.name}', app.width//2, 250, size=23)

def namePage_onKeyPress(app, key):
    if key=='enter':
        app.SMILES = smilesFromName(app.name)
        setActiveScreen('resultPage')
    elif key=='backspace':
        app.name = app.name[:-1]
    elif key == 'space':
        app.name += ' '
    else:
        app.name += key

def smilesFromName(name):
    url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}"+
    "/property/CanonicalSMILES/TXT")
    response = requests.get(url)
    
    if response.status_code == 200:
        smiles = response.text.strip()  
        L = smiles.splitlines()
        return L[0]
    else:
        return None


def commonPage_onScreenActivate(app):
    app.organicMolecules = [
    'Methane', 'Butane', 'Ethylene', 'Acetylene','Polythene', 'Polystyrene', 
    'Ethy-bromine', 'Chloroform', 'Methyl-alcohol', 'Firmament', 'Ethyl-alcohol' 
    ,'Glycerol', 'Formaldehyde', 'Acetaldehyde', 'Acetone', 'Formic-acid', 
    'Acetic-acid', 'Acetyl-chloride', 'Acetic-anhydride', 'Acetamide', 
    'Ethyl-acetate', 'Urea', 'Oxalic-acid', 'Glucose', 'Benzene', 'Toluene', 
    'Chloro-benzene', 'Nitro-benzene', 'Aniline', 'Phenol', 'Benzaldehyde', 
    'Benzoic acid', 'Ether', 'Carbon-tetrachloride', 'Urotropin', 'Gammexene']
    app.organicMolecules.sort()

    app.vitamins = ['Vitamin A', 'Vitamin B1', 'Vitamin B2', 'Vitamin B3', 
                    'Vitamin B5', 'Vitamin B6', 'Vitamin B7', 'Vitamin B9', 
                    'Vitamin B12', 'Vitamin C', 'Vitamin D', 'Vitamin E', 
                    'Vitamin E', 'Vitamin K']
    
    app.aminoAcids = [
    {"name": "Alanine", "threeLetter": "Ala"},
    {"name": "Arginine", "threeLetter": "Arg"},
    {"name": "Asparagine", "threeLetter": "Asn"},
    {"name": "Aspartic Acid", "threeLetter": "Asp"},
    {"name": "Cysteine", "threeLetter": "Cys"},
    {"name": "Glutamic Acid", "threeLetter": "Glu"},
    {"name": "Glutamine", "threeLetter": "Gln"},
    {"name": "Glycine", "threeLetter": "Gly"},
    {"name": "Histidine", "threeLetter": "His"},
    {"name": "Isoleucine", "threeLetter": "Ile"},
    {"name": "Leucine", "threeLetter": "Leu"},
    {"name": "Lysine", "threeLetter": "Lys"},
    {"name": "Methionine", "threeLetter": "Met"},
    {"name": "Phenylalanine", "threeLetter": "Phe"},
    {"name": "Proline", "threeLetter": "Pro"},
    {"name": "Serine", "threeLetter": "Ser"},
    {"name": "Threonine", "threeLetter": "Thr"},
    {"name": "Tryptophan", "threeLetter": "Trp"},
    {"name": "Tyrosine", "threeLetter": "Tyr"},
    {"name": "Valine", "threeLetter": "Val"}]

    app.aminoAcidsRow = 0
    app.vitaminRow = 0
    app.organicRow = 0


def commonPage_redrawAll(app):
    #amino acids
    drawLabel('Amino Acids', 150, 25, size=25, bold=True)
    aminoAcidsBoard = Board(10, 1, 50, 100, 200, 450)
    aminoAcidsBoard.drawBoard()
    drawRegularPolygon(150, 75, 18, 3, fill=None, border=rgb(58, 100, 136), 
                        borderWidth=3)
    drawRegularPolygon(150, 575, 18, 3, rotateAngle=180, fill=None, 
                        border=rgb(58, 100, 136), borderWidth=3)
    for row in range(10):
        drawLabel(f'{app.aminoAcids[row+app.aminoAcidsRow]['name']}', 150, 
                    122.5+45*row, size=20, bold=True)
    #vitamins
    drawLabel('Vitamins', 400, 25, size=25, bold=True)
    vitaminBoard = Board(10, 1, 300, 100, 200, 450)
    vitaminBoard.drawBoard()
    drawRegularPolygon(400, 75, 18, 3, fill=None, border=rgb(58, 100, 136), 
                        borderWidth=3)
    drawRegularPolygon(400, 575, 18, 3, rotateAngle=180, fill=None, 
                    border=rgb(58, 100, 136), borderWidth=3)
    for row in range(10):
        drawLabel(f'{app.vitamins[row+app.vitaminRow]}', 400, 122.5+45*row, 
                    size=20, bold=True)
    #organic molecules
    drawLabel('Ogranic Compounds', 650, 25, size=25, bold=True)
    organicMoleculesBoard = Board(10, 1, 550, 100, 200, 450)
    organicMoleculesBoard.drawBoard()
    drawRegularPolygon(650, 75, 18, 3, fill=None, border=rgb(58, 100, 136), 
                        borderWidth=3)
    drawRegularPolygon(650, 575, 18, 3, rotateAngle=180, fill=None, 
                    border=rgb(58, 100, 136), borderWidth=3)
    for row in range(10):
        drawLabel(f'{app.organicMolecules[row+app.organicRow]}', 650, 
                122.5+45*row, size=20, bold=True)
    
def commonPage_onMousePress(app, mouseX, mouseY):
    #click rows
    for y in range(100, 549, 45):
        if y<=mouseY<=y+45 and 50<=mouseX<=250:
            row = (y-100)//45
            app.SMILES = smilesFromName(app.aminoAcids[row+app.aminoAcidsRow]['name'])
            setActiveScreen('resultPage')
        elif y<=mouseY<=y+45 and 300<=mouseX<=500:
            row = (y-100)//45
            app.SMILES = smilesFromName(app.vitamins[row+app.vitaminRow])
            setActiveScreen('resultPage')
        elif y<=mouseY<=y+45 and 550<=mouseX<=750:
            row = (y-100)//45
            app.SMILES = smilesFromName(app.organicMolecules[row+app.organicRow])
            setActiveScreen('resultPage')
        
    #buttons
    if 137.5<=mouseX<=162.5:
        if 62.5<=mouseY<=87.5:
            app.aminoAcidsRow -= 1
            if app.aminoAcidsRow < 0:
                app.aminoAcidsRow = 0
        elif 562.5<=mouseY<=587.5:
            app.aminoAcidsRow += 1
            if app.aminoAcidsRow >= len(app.aminoAcids)-10:
                app.aminoAcidsRow = len(app.aminoAcids)-10
    elif 387.5<=mouseX<=412.5:
        if 62.5<=mouseY<=87.5:
            app.vitaminRow -= 1
            if app.vitaminRow < 0:
                app.vitaminRow = 0
        elif 562.5<=mouseY<=587.5:
            app.vitaminRow += 1
            if app.vitaminRow >= len(app.vitamins)-10:
                app.vitaminRow = len(app.vitamins)-10
    elif 637<=mouseX<=662.5:
        if 62.5<=mouseY<=87.5:
            app.organicRow -= 1
            if app.organicRow < 0:
                app.organicRow = 0
        elif 562.5<=mouseY<=587.5:
            app.organicRow += 1
            if app.organicRow >= len(app.organicMolecules)-10:
                app.organicRow = len(app.organicMolecules)-10

def generateSmiles(app):
    if len(app.selectedElements) == 0:
        return None
    startElement = app.selectedElements[0]  
    smiles = generateSmilesHelper(app, startElement)
    return smiles

def resultPage_onScreenActivate(app):
    if app.SMILES == None:
        app.SMILES = generateSmiles(app)
        app.stereo = smilesWithStereochemistry(app)
    app.stereoInfo = findIsomerInformation(app.SMILES)
    app.propertyInfo = getPropertyInfo(app.SMILES) 
    app.isomerNum = -1
    app.moleculeImage = generateImage(app.SMILES)
    app.isomerImages = generateIsomerImages(app.SMILES)
    app.presentImage = app.moleculeImage
    app.width = 800
    app.height = 600
    
def findIsomerInformation(smiles):
    mol = Chem.MolFromSmiles(smiles)

    chiralCenters = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    numChiralityCenters = len(chiralCenters)

    numDoubleBonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            numDoubleBonds += 1

    options = StereoEnumerationOptions(onlyUnassigned=False)
    isomers = list(EnumerateStereoisomers(mol, options=options))
    numStereoisomers = len(isomers)
    return [numDoubleBonds, numChiralityCenters, numStereoisomers]

def generateIsomerImages(smiles):
    mol = Chem.MolFromSmiles(smiles)
    options = StereoEnumerationOptions(onlyUnassigned=False)
    isomers = list(EnumerateStereoisomers(mol, options=options))

    results = []
    for i in range(len(isomers)):
        isomer = isomers[i]
        Draw.MolToFile(isomer, f"image{i+1}.png")
        results.append(f'image{i+1}.png')
    return results


def resultPage_redrawAll(app):
    #drawLine(0, 400, 800, 400)
    #drawRect(0, 0, 800, 800, fill=rgb(134, 176, 188))
    if app.propertyInfo[0] == None:
        drawLabel("Molecule Doesn't Exist", app.width//2, 25, size=35)
    else:
        drawLabel(f'IUPAC Name: {app.propertyInfo[0]}', app.width//2, 25, size=35)
        drawRect(25, 65, 750, 510, fill=rgb(201, 228, 231), border=rgb(38, 84, 124))
        drawLabel('Molecular Information', app.width//2, 70, size=28, 
                align='top')
        drawLabel(f'Molecular Formula: {app.propertyInfo[1]}', 35, 110, 
                align='left-top', size=25)
        drawLabel(f'Molecular Weight: {app.propertyInfo[2]} g/mol', 35, 140, 
                align='left-top', size=25)
        
        drawLabel('Hydrogen Bonding Information', app.width//2, 185, size=28, 
                align='top')
        drawLabel(f'H-Bond Donors: {app.propertyInfo[3]}', 35, 225, 
                align='left-top', size=25)
        drawLabel(f'H-Bond Acceptors: {app.propertyInfo[4]}', 35, 255, 
                align='left-top', size=25)

        drawLabel('Stereoisomer Information', app.width//4+30, 300, 
                align='top', size=28)
        drawLabel(f'Double bonds: {app.stereoInfo[0]}', 35, 340, 
                align='left-top', size=25)
        drawLabel(f'Chirality: {app.stereoInfo[1]}', 35, 370, 
                align='left-top', size=25)
        drawLabel(f'Total stereoisomers: {app.stereoInfo[2]}', 35, 400, 
                align='left-top', size=25)
        if app.stereoInfo[2]>1:
            drawRect(80, 450, 240, 50, fill=None, border=rgb(38, 84, 124))
            drawLabel("Click to see next isomer", 200, 475, size=20)
            drawLabel(f'{app.isomerNum+1}/{app.stereoInfo[2]} isomers', 
                    200, 510, size=18)
        if app.isomerNum==-1:
            drawImage(app.moleculeImage, 395, 300, width=360, height=240)
            drawRect(395, 300, 360, 240, fill=None, border=rgb(38, 84, 124))
        else:
            drawImage(app.isomerImages[app.isomerNum], 395, 300, width=360, height=240)
            drawRect(395, 300, 360, 240, fill=None, border=rgb(38, 84, 124))
        drawLabel('click to enlarge image', 575, 555, size=15)
    #home button image from https://www.flaticon.com/free-icon/home-button_61972
    drawImage('homeButton.png', 40, 510, width=50, height=50)
    drawRect(40, 510, 50, 50, fill=None, border=rgb(38, 84, 124))

def resultPage_onMousePress(app, mouseX, mouseY):
    if 395<=mouseX<=755 and 300<=mouseY<=540:
        setActiveScreen('imagePage')
    elif 80<=mouseX<=320 and 450<=mouseY<=500:
        app.isomerNum += 1
        if app.isomerNum == app.stereoInfo[2]:
            app.isomerNum = 0
        app.presentImage = app.isomerImages[app.isomerNum]
    elif 40<=mouseX<=90 and  510<=mouseY<=560:
        setActiveScreen('modeSelection')

def imagePage_onScreenActivate(app):
    imageWidth, imageHeight = getImageSize(app.presentImage)
    app.height = imageWidth*2
    app.width = imageHeight*2

def imagePage_redrawAll(app):
    drawImage(app.presentImage, 0, 0, height=app.height , width=app.width)
    drawLabel("click 'space' to return", app.width//2, app.height-20, size=15)

def imagePage_onKeyPress(app, key):
    setActiveScreen('resultPage')

def main():
    runAppWithScreens(initialScreen='start')

main() 