09-Backup
================
Kathleen Durkin
2025-06-09

- <a href="#01-archiving-a-repo-with-zenodo"
  id="toc-01-archiving-a-repo-with-zenodo">0.1 Archiving a repo with
  Zenodo</a>
  - <a href="#011-1-create-a-zenodo-account"
    id="toc-011-1-create-a-zenodo-account">0.1.1 1. Create a Zenodo
    account.</a>
  - <a href="#012-2-link-github-with-zenodo"
    id="toc-012-2-link-github-with-zenodo">0.1.2 2. Link GitHub with
    Zenodo.</a>
  - <a href="#013-3-authorize-zenodo-on-github"
    id="toc-013-3-authorize-zenodo-on-github">0.1.3 3. Authorize Zenodo on
    GitHub.</a>
  - <a href="#014-4-select-the-repository"
    id="toc-014-4-select-the-repository">0.1.4 4. Select the repository.</a>
  - <a href="#015-5-create-a-new-release-on-github"
    id="toc-015-5-create-a-new-release-on-github">0.1.5 5. Create a new
    release on GitHub.</a>
  - <a href="#016-6-check-zenodo-for-the-deposit"
    id="toc-016-6-check-zenodo-for-the-deposit">0.1.6 6. Check Zenodo for
    the deposit.</a>
  - <a href="#017-7-publish-the-archive-on-zenodo"
    id="toc-017-7-publish-the-archive-on-zenodo">0.1.7 7. Publish the
    archive on Zenodo.</a>

Full assignment details
[here](https://sr320.github.io/course-fish546-2025/assignments/09-backup.html).

Zenodo Instructions (provided
[here](https://sr320.github.io/course-fish546-2023/lectures/09-project.html#archiving-a-repo-with-zenodo)):

## 0.1 Archiving a repo with Zenodo

Archiving your GitHub repository with Zenodo allows you to have a
digital object identifier (DOI) for your repo, which is useful for
citing your software in academic papers. Here’s a step-by-step guide:

### 0.1.1 1. Create a Zenodo account.

You’ll need to visit Zenodo’s homepage (<https://zenodo.org/>) and sign
up for an account if you don’t already have one.

**Created an account under the username `shedurkin`, and linked it to my
ORCID ID.**

### 0.1.2 2. Link GitHub with Zenodo.

After you’ve set up your Zenodo account, you can navigate to the
“GitHub” section in the dashboard and link your GitHub account to
Zenodo.

### 0.1.3 3. Authorize Zenodo on GitHub.

Zenodo will request access to your GitHub repositories. You can either
allow Zenodo to access all your repositories or only select ones. After
you’ve made your choice, click “Authorize Zenodo.”

### 0.1.4 4. Select the repository.

Back in Zenodo, you should see a list of all your GitHub repositories.
You can choose to archive any repository by toggling the switch to “on”
next to the repository name. Zenodo will now create a “webhook” for this
repository, which will notify Zenodo whenever there is a new release of
the repository on GitHub.

### 0.1.5 5. Create a new release on GitHub.

If you haven’t done so already, you should create a new release of your
repository on GitHub. Go to your repository page, click “releases” then
“create a new release”. Tag version with the version number (for example
“v1.0”) and add some description about this release. Zenodo will be
notified once the release is published.

### 0.1.6 6. Check Zenodo for the deposit.

After creating a new release, go back to Zenodo. You should see the new
release in the “Upload” section. Zenodo automatically fills in some of
the metadata for you, such as the title and authors, but you can edit
this if needed.

### 0.1.7 7. Publish the archive on Zenodo.

Finally, you can publish the archive on Zenodo. Make sure all the
information is correct, and then click the “Publish” button. Your
repository is now archived, and you will have a DOI that you can use to
cite your software.

Remember that Zenodo will archive a new version of your repository every
time you create a new release on GitHub, so it’s important to create
meaningful and well-documented releases. Zenodo also allows you to link
different versions of your software together, so that it’s clear how the
software has changed over time.
