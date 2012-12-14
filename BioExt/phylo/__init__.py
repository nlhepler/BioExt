
from __future__ import division, print_function

from collections import deque, namedtuple


__all__ = ['Phylo']


class Clade(namedtuple('Clade', ['branch_length', 'children', 'depth', 'name', 'parent'])):

    def __init__(self, *args, **kwargs):
        super(Clade, self).__init__(*args, **kwargs)
        self.index = None

    def __repr__(self):
        branch_length = 1 if self.branch_length is None else self.branch_length
        attrs = ['branch_length=%g' % branch_length]
        if self.name:
            attrs.append("name='%s'" % self.name)
        return "Clade(%s)" % ', '.join(attrs)


class NewickParser(object):

    def __init__(self):
        pass

    @staticmethod
    def parse_branch(txt, children, parent, depth):
        try:
            if txt.find(':') < 0:
                name = txt
                bl = ''
            else:
                name, bl = txt.split(':')
            branch_length = float(bl) if len(bl) else None
            return Clade(branch_length, children, depth, name, parent)
        except ValueError:
            msg = "invalid Branch: '%s'" % txt
            raise RuntimeError(msg)

    @classmethod
    def parse_tree(cls, txt, parent, depth):
        BRANCH, BRANCHSET, ROOT = range(3)
        state = ROOT
        lst = deque()
        lwr = -1
        npar = 0
        childtxt = None
        for i, char in enumerate(txt):
            if state == ROOT:
                if char == '(':
                    state = BRANCHSET
                    npar = 1
                    lwr = i + 1
                elif char == ')':
                    raise RuntimeError("unbalanced paren: '%s'" % txt[:i+1])
                elif char == ';':
                    break
                else:
                    state = BRANCH
                    lwr = i + (1 if char == ',' else 0)
            elif state == BRANCHSET:
                if char == '(':
                    npar += 1
                if char == ')':
                    npar -= 1
                    if npar == 0:
                        childtxt = txt[lwr:i]
                        state = ROOT
            elif state == BRANCH:
                if char in ',;':
                    # we have to preinitialize the list so we can initialize
                    # the attributes in order to play by the read-only rules
                    # of namedtuple
                    children = []
                    node = cls.parse_branch(txt[lwr:i], children, parent, depth)
                    if childtxt is not None:
                        children.extend(cls.parse_tree(childtxt, node, depth + 1))
                        lst.extendleft(children[::-1])
                        childtxt = None
                    lst.appendleft(node)
                    state = ROOT
                elif char in '()':
                    raise RuntimeError("unbalanced paren: '%s'" % txt[lwr:i+1])

        # finalize
        if state == BRANCH:
            children = []
            node = cls.parse_branch(txt[lwr:], children, parent, depth)
            if childtxt is not None:
                children.extend(cls.parse_tree(childtxt, node, depth + 1))
                lst.extendleft(children[::-1])
                childtxt = None
            lst.appendleft(node)
        elif state == BRANCHSET:
            raise RuntimeError("unbalanced paren: '%s'" % txt[lwr:])

        # if the whole thing was in parenthesis
        if childtxt is not None:
            if len(lst):
                raise RuntimeError('something unexpected happened')
            return cls.parse_tree(childtxt, parent, depth)

        return lst

    @classmethod
    def parse(cls, txt):
        tree = cls.parse_tree(txt, None, 0)
        tree.reverse()
        return Tree(tree)


class Tree(deque):

    def __init__(self, *args, **kwargs):
        super(Tree, self).__init__(*args, **kwargs)
        cladogram = True
        for i, node in enumerate(self):
            node.index = i
            if node.branch_length is not None:
                cladogram = False
        self.cladogram = cladogram

    @staticmethod
    def __mrca(default, shared, node, acc=0):
        d = default if node.branch_length is None else node.branch_length
        if node.index in shared:
            return node.index, acc
        elif node.parent is None:
            return -1, acc + d
        else:
            return Tree.__mrca(default, shared, node.parent, acc + d)

    @staticmethod
    def __path_to_root(default, node, acc=0):
        d = default if node.branch_length is None else node.branch_length
        if node.parent is None:
            lst = {-1: acc + d}
        else:
            lst = Tree.__path_to_root(default, node.parent, acc + d)
        lst[node.index] = acc
        return lst

    def distance(self, node1, node2):
        d = 1 if self.cladogram else 0.0
        path1 = Tree.__path_to_root(d, node1)
        idx, path2 = Tree.__mrca(d, path1, node2)
        return path1[idx] + path2
        # shared = path1.keys() & path2.keys()
        # return min(path1[k] + path2[k] for k in shared)

    def get_terminals(self):
        return [n for n in self if len(n.children) == 0]

    def __repr__(self):
        lines = []
        if len(self):
            lines.extend((' ' * 4 * (node.depth + 1)) + repr(node) for node in self)
        lines.reverse()
        return 'Tree(weight=1.0 rooted=False)\n' + '\n'.join(lines)

    def __str__(self):
        return repr(self)


class Phylo(object):

    @staticmethod
    def read(fh, fmt):
        if fmt != 'newick':
            msg = 'only newick formatted trees supported at this time'
            raise ValueError(msg)
        return NewickParser.parse(fh.read())
